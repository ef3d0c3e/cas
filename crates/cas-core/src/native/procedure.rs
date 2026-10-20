//! Procedures (CAS-DOC-7 §8): `Block[{locals…}, statements…]`, run one
//! statement at a time in a frame that holds the values of its locals.
//!
//! Before a statement is evaluated, the values of the locals assigned so far
//! are put in place of their names, so nothing evaluated ever contains a
//! local: every part is a function of the state and of the node, and the
//! stamps and the memo stay sound with no special case (CAS-DOC-4 §2.3). The
//! block itself is held, and its result memoised as any native's. Lowering
//! checks that a local is assigned before it is read (O5), and that
//! parameters are not assigned: the rule put their values in the body.
//!
//! The statements:
//!
//! - `Set[x, e]`, `Set[{a, b}, e]`: `e` evaluated, then assigned to locals;
//! - `If[c, then, else]` whose branches are `Statements[…]`, groups run in
//!   the same frame; its condition must be `True` or `False`, as in Xcas
//!   ([`EvalError::UndecidedCondition`]);
//! - `For[k, r, body]`: `r` a `Range[a, b]` of integers (empty when `b < a`),
//!   a list or a set ([`EvalError::NotIterable`] otherwise), `k` taking each
//!   value in turn; `While[c, body]`;
//! - `Return[e]`, `Return[]`, `Break[]`, `Continue[]`;
//! - anything else, an expression, evaluated.
//!
//! The value of a group is that of its last statement: an expression's
//! value, the value assigned, a branch's value, `Null` for a loop. A loop
//! checks the interrupt at each round, since it may run without a rewrite
//! (`while true do end`). The block and its nested groups take levels of
//! `Limits::recursion`, so the frames of deep procedures stay within the
//! stack.

use std::sync::atomic::Ordering;

use foldhash::fast::FixedState;
use hashbrown::{HashMap, HashSet};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// How a statement ends.
enum Flow {
    /// On to the next statement, with its value.
    Next(Expr),
    /// Out of the block, with its value.
    Return(Expr),
    /// Out of the innermost loop.
    Break,
    /// To the next round of the innermost loop.
    Continue,
}

/// The values of the locals assigned so far.
struct Frame {
    locals: HashSet<Symbol, FixedState>,
    values: HashMap<Symbol, Expr, FixedState>,
}

fn has_head(e: &Expr, s: Symbol) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(s)
}

/// `Block[{locals…}, statements…]`.
pub(super) fn block(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [locals, statements @ ..] = e.args() else {
        return Ok(None);
    };
    if !has_head(locals, Symbol::LIST) {
        return Ok(None);
    }
    let Some(locals) = locals.args().iter().map(Expr::as_symbol).collect() else {
        return Ok(None);
    };
    let mut frame = Frame {
        locals,
        values: HashMap::with_hasher(FixedState::default()),
    };
    // The frames of the interpreter take a level, as the matcher's do.
    ev.enter()?;
    let flow = run(ev, &mut frame, statements);
    ev.leave();
    Ok(Some(match flow? {
        Flow::Next(value) | Flow::Return(value) => value,
        Flow::Break | Flow::Continue => Symbol::NULL.into(),
    }))
}

/// Runs `statements` in order.
fn run(ev: &mut Eval<'_>, frame: &mut Frame, statements: &[Expr]) -> Result<Flow, EvalError> {
    let mut value = Expr::from(Symbol::NULL);
    for s in statements {
        match statement(ev, frame, s)? {
            Flow::Next(v) => value = v,
            other => return Ok(other),
        }
    }
    Ok(Flow::Next(value))
}

/// Runs a group of statements a level down.
fn group(ev: &mut Eval<'_>, frame: &mut Frame, body: &Expr) -> Result<Flow, EvalError> {
    ev.enter()?;
    let flow = if has_head(body, Symbol::STATEMENTS) {
        run(ev, frame, body.args())
    } else {
        statement(ev, frame, body)
    };
    ev.leave();
    flow
}

fn statement(ev: &mut Eval<'_>, frame: &mut Frame, s: &Expr) -> Result<Flow, EvalError> {
    let head = s.head().and_then(Expr::as_symbol);
    Ok(match (head, s.args()) {
        (Some(Symbol::SET), [target, value]) => {
            let value = ev.eval(&frame.substitute(value))?;
            frame.assign(target, &value)?;
            Flow::Next(value)
        }
        (Some(Symbol::IF), [c, then]) => {
            if decide(ev, frame, c)? {
                group(ev, frame, then)?
            } else {
                Flow::Next(Symbol::NULL.into())
            }
        }
        (Some(Symbol::IF), [c, then, other]) => {
            let branch = if decide(ev, frame, c)? { then } else { other };
            group(ev, frame, branch)?
        }
        (Some(Symbol::FOR), [variable, values, body]) => {
            let Some(variable) = variable.as_symbol() else {
                return Err(EvalError::NotIterable);
            };
            let values = ev.eval(&frame.substitute(values))?;
            for_each(ev, frame, variable, &values, body)?
        }
        (Some(Symbol::WHILE), [c, body]) => loop {
            interrupted(ev)?;
            if !decide(ev, frame, c)? {
                break Flow::Next(Symbol::NULL.into());
            }
            match group(ev, frame, body)? {
                Flow::Next(_) | Flow::Continue => {}
                Flow::Break => break Flow::Next(Symbol::NULL.into()),
                Flow::Return(v) => break Flow::Return(v),
            }
        },
        (Some(Symbol::RETURN), []) => Flow::Return(Symbol::NULL.into()),
        (Some(Symbol::RETURN), [value]) => Flow::Return(ev.eval(&frame.substitute(value))?),
        (Some(Symbol::BREAK), []) => Flow::Break,
        (Some(Symbol::CONTINUE), []) => Flow::Continue,
        (Some(Symbol::STATEMENTS), _) => group(ev, frame, s)?,
        _ => Flow::Next(ev.eval(&frame.substitute(s))?),
    })
}

/// Runs `body` for each of `values` given to `variable`.
fn for_each(
    ev: &mut Eval<'_>,
    frame: &mut Frame,
    variable: Symbol,
    values: &Expr,
    body: &Expr,
) -> Result<Flow, EvalError> {
    let round = |ev: &mut Eval<'_>, frame: &mut Frame, value: Expr| {
        interrupted(ev)?;
        frame.values.insert(variable, value);
        group(ev, frame, body)
    };
    let done = Flow::Next(Symbol::NULL.into());
    if has_head(values, Symbol::RANGE) {
        let [from, to] = values.args() else {
            return Err(EvalError::NotIterable);
        };
        let (Some(from), Some(to)) = (from.to_integer(), to.to_integer()) else {
            return Err(EvalError::NotIterable);
        };
        let one = cas_num::Integer::from(1);
        let mut value = from;
        while value <= to {
            match round(ev, frame, Expr::from(&value))? {
                Flow::Next(_) | Flow::Continue => {}
                Flow::Break => break,
                Flow::Return(v) => return Ok(Flow::Return(v)),
            }
            value = &value + &one;
        }
        return Ok(done);
    }
    if !(has_head(values, Symbol::LIST) || has_head(values, Symbol::FINITE_SET)) {
        return Err(EvalError::NotIterable);
    }
    for value in values.args() {
        match round(ev, frame, value.clone())? {
            Flow::Next(_) | Flow::Continue => {}
            Flow::Break => break,
            Flow::Return(v) => return Ok(Flow::Return(v)),
        }
    }
    Ok(done)
}

/// Whether the condition `c` holds: it must evaluate to `True` or `False`.
fn decide(ev: &mut Eval<'_>, frame: &Frame, c: &Expr) -> Result<bool, EvalError> {
    match ev.eval(&frame.substitute(c))?.as_symbol() {
        Some(Symbol::TRUE) => Ok(true),
        Some(Symbol::FALSE) => Ok(false),
        _ => Err(EvalError::UndecidedCondition),
    }
}

/// [`EvalError::Interrupted`] once the evaluation was interrupted.
fn interrupted(ev: &Eval<'_>) -> Result<(), EvalError> {
    if ev.context().interrupt_flag().load(Ordering::Relaxed) {
        Err(EvalError::Interrupted)
    } else {
        Ok(())
    }
}

impl Frame {
    /// Assigns `value` to the local `target`, or to the locals of a list
    /// element by element.
    fn assign(&mut self, target: &Expr, value: &Expr) -> Result<(), EvalError> {
        match target.view() {
            View::Symbol(s) if self.locals.contains(&s) => {
                self.values.insert(s, value.clone());
                Ok(())
            }
            View::Apply { args, .. } if has_head(target, Symbol::LIST) => {
                if !has_head(value, Symbol::LIST) || value.args().len() != args.len() {
                    return Err(EvalError::Shape);
                }
                for (t, v) in args.iter().zip(value.args()) {
                    self.assign(t, v)?;
                }
                Ok(())
            }
            // Lowering makes every name assigned in a block one of its
            // locals; anything else is not written.
            _ => Err(EvalError::Shape),
        }
    }

    /// `e` with the values of the locals assigned so far in place of their
    /// names.
    fn substitute(&self, e: &Expr) -> Expr {
        if self.values.is_empty() {
            return e.clone();
        }
        substitute(e, &self.values)
    }
}

/// `e` with `values` in place of their names, except where a name is bound
/// again: by a nested block's locals, a function's parameters, or a pattern
/// variable, in its rule.
fn substitute(e: &Expr, values: &HashMap<Symbol, Expr, FixedState>) -> Expr {
    e.replace(|part| match part.view() {
        View::Symbol(s) => values.get(&s).cloned(),
        View::Apply { head, args } => {
            let bound = bound(head.as_symbol()?, args);
            if bound.is_empty() {
                return None;
            }
            let mut inner = values.clone();
            for s in &bound {
                inner.remove(s);
            }
            // The names bound stay, where they are declared too.
            Some(part.map_args(|a| substitute(a, &inner)))
        }
        _ => None,
    })
}

/// The names an application binds for its arguments.
fn bound(head: Symbol, args: &[Expr]) -> Vec<Symbol> {
    let names = |e: &Expr| -> Vec<Symbol> {
        match e.view() {
            View::Symbol(s) => vec![s],
            _ if has_head(e, Symbol::LIST) => e.args().iter().filter_map(Expr::as_symbol).collect(),
            _ => Vec::new(),
        }
    };
    match (head, args) {
        (Symbol::BLOCK, [locals, ..]) | (Symbol::FUNCTION, [locals, _]) => names(locals),
        (Symbol::PATTERN, [name, _]) => names(name),
        (Symbol::RULE | Symbol::RULE_DELAYED, [lhs, _]) => {
            let mut variables = Vec::new();
            pattern_variables(lhs, &mut variables);
            variables
        }
        _ => Vec::new(),
    }
}

/// The names `Pattern[name, …]` binds in `e`.
fn pattern_variables(e: &Expr, out: &mut Vec<Symbol>) {
    if let (Some(Symbol::PATTERN), [name, _]) = (e.head().and_then(Expr::as_symbol), e.args()) {
        out.extend(name.as_symbol());
    }
    for a in e.args() {
        pattern_variables(a, out);
    }
}

#[cfg(test)]
mod tests;
