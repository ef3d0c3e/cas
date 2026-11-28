//! A static index of the names in one text (CAS-DOC-60 §4): what each
//! [`Name`](crate::SyntaxKind::Name) node stands for, where it is defined and used,
//! and in which scope. It mirrors lowering's rules (R2.2-004) with a walk of
//! its own over the tree, without lowering, so it works on half-typed text:
//! a statement with a syntax error keeps the names its tree shows, and a
//! name inside an error node is none.
//!
//! A scope binds its keys before its contents are walked (a rule's
//! parameters, a lambda's, a block's locals as lowering finds them), so a
//! use resolves to the innermost scope that binds its key wherever it stands
//! in that scope. A reserved key always resolves to its target; any other
//! key to the one global of that key in the text.

use std::collections::HashMap;
use std::ops::Range;

use crate::analysis::cursor::token_at;
use crate::lower::{arguments, assigned_names, has_statements, name_of, range};
use crate::names::{Names, SyntaxForm, Target, TypeName};
use crate::parser::Parse;
use crate::syntax::{SyntaxKind as K, SyntaxNode};

/// The index of a symbol in [`Symbols::symbols`], valid for one [`Symbols`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SymbolId(pub usize);

/// The index of a scope in [`Symbols::scopes`], valid for one [`Symbols`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ScopeId(pub usize);

/// What one symbol is, and where it is defined and used.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SymbolInfo {
    /// A user symbol's key, as [`name`](crate::name) gives it (`pi` for
    /// `π`), unescaped; a reserved target's or a type's spelling.
    pub name: String,
    /// What it stands for.
    pub kind: SymbolKind,
    /// The scope that binds it: the document for globals, reserved names
    /// and types.
    pub scope: ScopeId,
    /// Its binding sites, in text order.
    pub definitions: Vec<Definition>,
    /// Its other occurrences, in text order.
    pub uses: Vec<Range<usize>>,
}

/// A binding site.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Definition {
    /// The name's bytes.
    pub name: Range<usize>,
    /// The construct's bytes: the `Definition` statement, the `LocalItem`,
    /// the `ForStmt`, the parameter, or the `Fill` call.
    pub node: Range<usize>,
}

/// A `Name` node that stands for a symbol.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Occurrence {
    /// The name's bytes.
    pub range: Range<usize>,
    /// What it stands for.
    pub symbol: SymbolId,
    /// Whether it is a binding site.
    pub definition: bool,
}

/// A region that binds names.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Scope {
    /// What kind of region.
    pub kind: ScopeKind,
    /// Its bytes.
    pub range: Range<usize>,
    /// The scope around it; `None` for the document.
    pub parent: Option<ScopeId>,
}

/// The kinds of [`Scope`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ScopeKind {
    /// The whole text, which holds the globals.
    Document,
    /// A definition `f(params) := body [when guard]`, binding its
    /// parameters.
    Rule,
    /// A `Lambda`, binding its parameters.
    Lambda,
    /// A `begin … end`, or an `if` with statements that is a definition's
    /// body, binding its locals.
    Block,
}

/// What a symbol stands for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SymbolKind {
    /// A builtin of cas-core, settings such as `Digits` included.
    Builtin(cas_core::Symbol),
    /// The imaginary unit `i`.
    Imaginary,
    /// A builtin that is reserved but not defined yet.
    Planned,
    /// A form that lowering or the session handles itself.
    Syntax(SyntaxForm),
    /// A type after `::`.
    Type(TypeName),
    /// A global with a rule definition.
    Function,
    /// A global without one.
    Variable,
    /// A parameter of a rule or a lambda.
    Parameter,
    /// A local of a block.
    Local,
}

impl SymbolKind {
    /// Whether it is a user's symbol: a function, a variable, a parameter
    /// or a local.
    #[must_use]
    pub fn is_user(self) -> bool {
        matches!(
            self,
            SymbolKind::Function | SymbolKind::Variable | SymbolKind::Parameter | SymbolKind::Local
        )
    }

    fn of(target: &Target) -> SymbolKind {
        match target {
            Target::Builtin(s) => SymbolKind::Builtin(*s),
            Target::Planned(_) => SymbolKind::Planned,
            Target::Imaginary => SymbolKind::Imaginary,
            Target::Syntax(form) => SymbolKind::Syntax(*form),
        }
    }
}

/// The names of one text (CAS-DOC-60 §4).
#[derive(Clone, Debug)]
pub struct Symbols {
    parse: Parse,
    infos: Vec<SymbolInfo>,
    occurrences: Vec<Occurrence>,
    scopes: Vec<Scope>,
}

impl Symbols {
    /// The symbols, indexed by [`SymbolId`].
    #[must_use]
    pub fn symbols(&self) -> &[SymbolInfo] {
        &self.infos
    }

    /// The occurrences, in text order and disjoint.
    #[must_use]
    pub fn occurrences(&self) -> &[Occurrence] {
        &self.occurrences
    }

    /// The scopes, indexed by [`ScopeId`]: the document first, each parent
    /// before its children.
    #[must_use]
    pub fn scopes(&self) -> &[Scope] {
        &self.scopes
    }

    /// The occurrence covering the token at `offset`, as
    /// [`token_at`] picks it.
    #[must_use]
    pub fn at(&self, offset: usize) -> Option<&Occurrence> {
        let token = token_at(&self.parse, offset)?;
        let r = token.text_range();
        let (start, end) = (usize::from(r.start()), usize::from(r.end()));
        let after = self.occurrences.partition_point(|o| o.range.start <= start);
        let o = &self.occurrences[after.checked_sub(1)?];
        (end <= o.range.end).then_some(o)
    }
}

/// The names of the text `parse` holds, with `names` for what is reserved.
#[must_use]
pub fn symbols(parse: &Parse, names: &Names) -> Symbols {
    let root = parse.syntax();
    let mut index = Index {
        names,
        infos: Vec::new(),
        occurrences: Vec::new(),
        scopes: Vec::new(),
        stack: Vec::new(),
        parameters: Vec::new(),
        globals: HashMap::new(),
        reserved: HashMap::new(),
        types: HashMap::new(),
    };
    index.open(ScopeKind::Document, range(&root), &[], SymbolKind::Variable);
    index.walk(&root);
    let mut occurrences = index.occurrences;
    occurrences.sort_by_key(|o| o.range.start);
    Symbols {
        parse: parse.clone(),
        infos: index.infos,
        occurrences,
        scopes: index.scopes,
    }
}

/// A name of a rule's left-hand side, classified before the rule's scope
/// opens.
enum Part {
    /// The user function the rule defines.
    Head(SyntaxNode),
    /// A parameter, with its construct's bytes.
    Parameter(SyntaxNode, Range<usize>),
    /// A name that is used.
    Use(SyntaxNode),
    /// A type after `::`.
    Type(SyntaxNode),
    /// A node walked as usual, inside the rule's scope.
    Walk(SyntaxNode),
}

struct Index<'a> {
    names: &'a Names,
    infos: Vec<SymbolInfo>,
    occurrences: Vec<Occurrence>,
    scopes: Vec<Scope>,
    /// The scopes around the walk, innermost last, with the keys they bind.
    stack: Vec<(ScopeId, HashMap<String, SymbolId>)>,
    /// Lowering's `parameters`: those of the user rules and lambdas around
    /// the walk, which a block's locals leave out.
    parameters: Vec<String>,
    globals: HashMap<String, SymbolId>,
    reserved: HashMap<Target, SymbolId>,
    types: HashMap<TypeName, SymbolId>,
}

impl Index<'_> {
    fn walk(&mut self, node: &SyntaxNode) {
        match node.kind() {
            K::ErrorNode => {}
            K::Name => self.name(node, None),
            K::Definition => self.definition(node),
            K::LocalItem | K::ForStmt => self.bound(node),
            K::Lambda => self.lambda(node),
            K::BlockExpr => {
                let statements: Vec<SyntaxNode> = node
                    .children()
                    .filter(|c| c.kind() == K::Body)
                    .flat_map(|body| body.children())
                    .collect();
                self.block(range(node), &statements);
                self.children(node);
                self.close();
            }
            K::ExprStmt => {
                if !node.first_child().is_some_and(|e| self.fill(&e)) {
                    self.children(node);
                }
            }
            K::Typed => {
                for (i, c) in node.children().enumerate() {
                    if i == 1 && c.kind() == K::Name {
                        self.type_name(&c);
                    } else {
                        self.walk(&c);
                    }
                }
            }
            _ => self.children(node),
        }
    }

    fn children(&mut self, node: &SyntaxNode) {
        for c in node.children() {
            self.walk(&c);
        }
    }

    /// A `LocalItem` or a `ForStmt`: its first child, a name, is a binding
    /// site of the whole node.
    fn bound(&mut self, node: &SyntaxNode) {
        for (i, c) in node.children().enumerate() {
            if i == 0 && c.kind() == K::Name {
                self.name(&c, Some(range(node)));
            } else {
                self.walk(&c);
            }
        }
    }

    /// `lhs := rhs [when guard]`: an assignment to a name, a list of names
    /// or an element, else a rule.
    fn definition(&mut self, node: &SyntaxNode) {
        let site = range(node);
        let colon = node
            .children_with_tokens()
            .find(|e| e.kind() == K::ColonEq)
            .map_or(site.end, |e| usize::from(e.text_range().start()));
        let (lhs, rest): (Vec<SyntaxNode>, Vec<SyntaxNode>) =
            node.children().partition(|c| range(c).end <= colon);
        let Some(lhs) = lhs.into_iter().next() else {
            return self.children(node);
        };
        match lhs.kind() {
            K::Name => self.name(&lhs, Some(site)),
            K::ListExpr => {
                for item in lhs.children() {
                    if item.kind() == K::Name {
                        self.name(&item, Some(site.clone()));
                    } else {
                        self.walk(&item);
                    }
                }
            }
            K::IndexExpr => {
                for (i, c) in lhs.children().enumerate() {
                    if i == 0 && c.kind() == K::Name {
                        self.name(&c, Some(site.clone()));
                    } else {
                        self.walk(&c);
                    }
                }
            }
            _ => return self.rule(&site, &lhs, &rest),
        }
        for c in &rest {
            self.walk(c);
        }
    }

    /// `f(params) := body`, or a definition of a builtin (an up-rule), whose
    /// `::`-tagged names and spreads are its parameters.
    fn rule(&mut self, site: &Range<usize>, lhs: &SyntaxNode, rest: &[SyntaxNode]) {
        let mut parts = Vec::new();
        let user = self.user_head(lhs);
        if user {
            self.defined_application(lhs, &mut parts);
        } else {
            self.literal(lhs, &mut parts);
        }
        let mut keys = Vec::new();
        for part in &parts {
            if let Part::Parameter(n, _) = part {
                push_new(&mut keys, name_of(n));
            }
        }
        self.open(ScopeKind::Rule, site.clone(), &keys, SymbolKind::Parameter);
        // An up-rule leaves lowering's `parameters` as they are.
        let saved = user.then(|| std::mem::replace(&mut self.parameters, keys));
        for part in parts {
            match part {
                Part::Head(n) => {
                    let id = self.global(name_of(&n));
                    self.infos[id.0].kind = SymbolKind::Function;
                    self.add(&n, id, Some(site.clone()));
                }
                Part::Parameter(n, node) => self.name(&n, Some(node)),
                Part::Use(n) => self.name(&n, None),
                Part::Type(n) => self.type_name(&n),
                Part::Walk(n) => self.walk(&n),
            }
        }
        for (i, c) in rest.iter().enumerate() {
            // An `if` with statements as the body is a block of its own.
            if i == 0 && c.kind() == K::IfExpr && has_statements(c) {
                self.block(range(c), std::slice::from_ref(c));
                self.children(c);
                self.close();
            } else {
                self.walk(c);
            }
        }
        if let Some(saved) = saved {
            self.parameters = saved;
        }
        self.close();
    }

    /// Whether `lhs` is an application whose innermost callee is a user
    /// name.
    fn user_head(&self, lhs: &SyntaxNode) -> bool {
        if lhs.kind() != K::Application {
            return false;
        }
        let mut callee = lhs.first_child();
        while let Some(c) = callee.as_ref().filter(|c| c.kind() == K::Application) {
            callee = c.first_child();
        }
        callee.is_some_and(|c| c.kind() == K::Name && !self.names.is_reserved(&name_of(&c)))
    }

    /// A user rule's left-hand side, `f(a)(b)`: the head, then the arguments
    /// of each level as parameters.
    fn defined_application(&self, node: &SyntaxNode, parts: &mut Vec<Part>) {
        for (i, c) in node.children().enumerate() {
            match c.kind() {
                K::Application if i == 0 => self.defined_application(&c, parts),
                K::Name if i == 0 => parts.push(Part::Head(c)),
                K::ArgList => self.parameters_in(&c, parts),
                _ => self.pattern(&c, range(&c), parts),
            }
        }
    }

    /// The arguments of `args`, each a parameter: `n = 1` binds `n`, its
    /// default is a value.
    fn parameters_in(&self, args: &SyntaxNode, parts: &mut Vec<Part>) {
        for arg in args.children() {
            let site = range(&arg);
            if arg.kind() == K::Equation {
                for (i, side) in arg.children().enumerate() {
                    if i == 0 {
                        self.pattern(&side, site.clone(), parts);
                    } else {
                        parts.push(Part::Walk(side));
                    }
                }
            } else {
                self.pattern(&arg, site, parts);
            }
        }
    }

    /// A pattern of a user rule, as lowering's `Parameters` mode reads it:
    /// every name not reserved is a parameter, but a callee is a value.
    fn pattern(&self, node: &SyntaxNode, site: Range<usize>, parts: &mut Vec<Part>) {
        match node.kind() {
            K::ErrorNode => {}
            K::Name if self.names.is_reserved(&name_of(node)) => {
                parts.push(Part::Use(node.clone()));
            }
            K::Name => parts.push(Part::Parameter(node.clone(), site)),
            K::Lambda | K::BlockExpr => parts.push(Part::Walk(node.clone())),
            _ => {
                for (i, c) in node.children().enumerate() {
                    match c.kind() {
                        K::Name if i == 1 && node.kind() == K::Typed => parts.push(Part::Type(c)),
                        K::Name if i == 0 && node.kind() == K::Application => {
                            parts.push(Part::Use(c));
                        }
                        K::ArgList => self.parameters_in(&c, parts),
                        _ => self.pattern(&c, site.clone(), parts),
                    }
                }
            }
        }
    }

    /// The left-hand side of an up-rule, literal: a `::`-tagged name or a
    /// spread is a parameter, any other name a use.
    fn literal(&self, node: &SyntaxNode, parts: &mut Vec<Part>) {
        match node.kind() {
            K::ErrorNode => {}
            K::Name => parts.push(Part::Use(node.clone())),
            K::Lambda | K::BlockExpr => parts.push(Part::Walk(node.clone())),
            kind => {
                let tags = matches!(kind, K::Typed | K::Spread);
                for (i, c) in node.children().enumerate() {
                    match c.kind() {
                        K::Name if i == 0 && tags && !self.names.is_reserved(&name_of(&c)) => {
                            parts.push(Part::Parameter(c, range(node)));
                        }
                        K::Name if i == 1 && kind == K::Typed => parts.push(Part::Type(c)),
                        _ => self.literal(&c, parts),
                    }
                }
            }
        }
    }

    /// `x -> body`, `(x, y) -> body`.
    fn lambda(&mut self, node: &SyntaxNode) {
        let mut children = node.children();
        let Some(params) = children.next() else {
            return;
        };
        let names: Vec<SyntaxNode> = if params.kind() == K::Name {
            vec![params]
        } else {
            params.children().filter(|n| n.kind() == K::Name).collect()
        };
        let mut keys = Vec::new();
        for n in &names {
            let key = name_of(n);
            if !self.names.is_reserved(&key) {
                push_new(&mut keys, key);
            }
        }
        self.open(ScopeKind::Lambda, range(node), &keys, SymbolKind::Parameter);
        let outer = self.parameters.len();
        self.parameters.extend(keys);
        for n in &names {
            self.name(n, Some(range(n)));
        }
        for c in children {
            self.walk(&c);
        }
        self.parameters.truncate(outer);
        self.close();
    }

    /// Opens a block of `statements`, binding lowering's locals: the names
    /// they assign, minus reserved names and the parameters around.
    fn block(&mut self, range: Range<usize>, statements: &[SyntaxNode]) {
        let mut locals = Vec::new();
        for s in statements {
            assigned_names(s, &mut locals);
        }
        locals.retain(|n| !self.names.is_reserved(n) && !self.parameters.contains(n));
        self.open(ScopeKind::Block, range, &locals, SymbolKind::Local);
    }

    /// `Fill(v, A)`, which assigns `A`: whether `node` is one.
    fn fill(&mut self, node: &SyntaxNode) -> bool {
        if node.kind() != K::Application {
            return false;
        }
        let Some(callee) = node.first_child().filter(|c| c.kind() == K::Name) else {
            return false;
        };
        if self.names.target(&name_of(&callee)) != Some(&Target::Syntax(SyntaxForm::Fill)) {
            return false;
        }
        let args = arguments(node);
        let [v, a] = &args[..] else {
            return false;
        };
        if a.kind() != K::Name {
            return false;
        }
        self.name(&callee, None);
        self.walk(v);
        self.name(a, Some(range(node)));
        true
    }

    /// An occurrence of the name `node`, a binding site of `site` unless it
    /// is reserved.
    fn name(&mut self, node: &SyntaxNode, site: Option<Range<usize>>) {
        let key = name_of(node);
        if let Some(target) = self.names.target(&key) {
            let id = if let Some(&id) = self.reserved.get(target) {
                id
            } else {
                let name = self.names.spelling(target).unwrap_or(&key).to_owned();
                let id = self.symbol(name, SymbolKind::of(target), ScopeId(0));
                self.reserved.insert(target.clone(), id);
                id
            };
            return self.add(node, id, None);
        }
        let bound = self
            .stack
            .iter()
            .rev()
            .find_map(|(_, keys)| keys.get(&key).copied());
        let id = bound.unwrap_or_else(|| self.global(key));
        self.add(node, id, site);
    }

    /// A type name after `::`: an occurrence when the table knows it.
    fn type_name(&mut self, node: &SyntaxNode) {
        let key = name_of(node);
        let Some(ty) = self.names.type_name(&key) else {
            return;
        };
        let id = if let Some(&id) = self.types.get(&ty) {
            id
        } else {
            let name = self.names.type_spelling(ty).unwrap_or(&key).to_owned();
            let id = self.symbol(name, SymbolKind::Type(ty), ScopeId(0));
            self.types.insert(ty, id);
            id
        };
        self.add(node, id, None);
    }

    /// The global of `key`.
    fn global(&mut self, key: String) -> SymbolId {
        if let Some(&id) = self.globals.get(&key) {
            return id;
        }
        let id = self.symbol(key.clone(), SymbolKind::Variable, ScopeId(0));
        self.globals.insert(key, id);
        id
    }

    fn symbol(&mut self, name: String, kind: SymbolKind, scope: ScopeId) -> SymbolId {
        self.infos.push(SymbolInfo {
            name,
            kind,
            scope,
            definitions: Vec::new(),
            uses: Vec::new(),
        });
        SymbolId(self.infos.len() - 1)
    }

    fn add(&mut self, node: &SyntaxNode, id: SymbolId, site: Option<Range<usize>>) {
        let name = range(node);
        let definition = site.is_some();
        let info = &mut self.infos[id.0];
        match site {
            Some(node) => info.definitions.push(Definition {
                name: name.clone(),
                node,
            }),
            None => info.uses.push(name.clone()),
        }
        self.occurrences.push(Occurrence {
            range: name,
            symbol: id,
            definition,
        });
    }

    /// Opens a scope binding `keys`, each a new symbol of `kind`.
    fn open(&mut self, kind: ScopeKind, range: Range<usize>, keys: &[String], symbol: SymbolKind) {
        let id = ScopeId(self.scopes.len());
        self.scopes.push(Scope {
            kind,
            range,
            parent: self.stack.last().map(|(id, _)| *id),
        });
        let mut bound = HashMap::new();
        for key in keys {
            if !bound.contains_key(key) {
                let s = self.symbol(key.clone(), symbol, id);
                bound.insert(key.clone(), s);
            }
        }
        self.stack.push((id, bound));
    }

    fn close(&mut self) {
        self.stack.pop();
    }
}

/// Pushes `key` unless `keys` has it.
fn push_new(keys: &mut Vec<String>, key: String) {
    if !keys.contains(&key) {
        keys.push(key);
    }
}

#[cfg(test)]
mod props;
#[cfg(test)]
mod tests;
