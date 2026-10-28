//! Token and node kinds, and the glue to `rowan`'s lossless trees.

/// Declares [`SyntaxKind`] and the table that turns rowan's raw kinds back
/// into it.
macro_rules! syntax_kinds {
    ($($(#[$doc:meta])* $kind:ident,)*) => {
        /// The kind of a token or of a node of the syntax tree.
        ///
        /// Aliases share a kind with what they stand for: `×` is a [`Star`]
        /// and `≤` a [`Le`], with their own text (CAS-DOC-7 §1). `π`, `ℯ` and
        /// `∞` are [`Ident`]s of one character; [`name`](crate::name) gives
        /// `pi`, `e` and `inf` for them.
        ///
        /// [`Star`]: SyntaxKind::Star
        /// [`Le`]: SyntaxKind::Le
        /// [`Ident`]: SyntaxKind::Ident
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
        #[repr(u16)]
        pub enum SyntaxKind {
            $($(#[$doc])* $kind,)*
        }

        impl SyntaxKind {
            /// Every kind, indexed by its raw value.
            const ALL: &[SyntaxKind] = &[$(SyntaxKind::$kind,)*];
        }
    };
}

syntax_kinds! {
    /// Spaces, tabs and any other whitespace except line breaks.
    Whitespace,
    /// One line break: `\n`, `\r\n` or `\r`. Whether it ends a statement is
    /// the parser's decision (CAS-DOC-7 §4).
    Newline,
    /// `// …`, up to the line break.
    LineComment,
    /// `/* … */`, which nests.
    BlockComment,

    /// A number: `42`, `1_000`, `0x1f`, `0.1`, `1.5e-3`.
    Number,
    /// A string, quotes included.
    String,
    /// An exponent in superscript digits, with an optional sign: `²`, `⁻¹`.
    Superscript,
    /// A name, or one of the standalone `π`, `ℯ`, `∞`.
    Ident,

    /// `if`.
    IfKw,
    /// `then`.
    ThenKw,
    /// `elif`.
    ElifKw,
    /// `else`.
    ElseKw,
    /// `end`.
    EndKw,
    /// `for`.
    ForKw,
    /// `in`, `∈`.
    InKw,
    /// `do`.
    DoKw,
    /// `while`.
    WhileKw,
    /// `begin`.
    BeginKw,
    /// `return`.
    ReturnKw,
    /// `break`.
    BreakKw,
    /// `continue`.
    ContinueKw,
    /// `local`.
    LocalKw,
    /// `when`.
    WhenKw,
    /// `and`.
    AndKw,
    /// `or`.
    OrKw,
    /// `not`.
    NotKw,
    /// `union`, `∪` (CAS-DOC-32 §3).
    UnionKw,
    /// `intersect`, `∩`.
    IntersectKw,
    /// `minus`, `∖`.
    MinusKw,

    /// `+`.
    Plus,
    /// `-`, `−`.
    Minus,
    /// `*`, `×`, `·`, `⋅`.
    Star,
    /// `/`, `÷`.
    Slash,
    /// `^`.
    Caret,
    /// `.*`.
    DotStar,
    /// `./`.
    DotSlash,
    /// `.^`.
    DotCaret,
    /// `!`.
    Bang,
    /// `=`.
    Eq,
    /// `==`.
    EqEq,
    /// `!=`, `≠`.
    Neq,
    /// `<`.
    Lt,
    /// `<=`, `≤`.
    Le,
    /// `>`.
    Gt,
    /// `>=`, `≥`.
    Ge,
    /// `..`.
    DotDot,
    /// `...`.
    Ellipsis,
    /// `::`.
    ColonColon,
    /// `:=`.
    ColonEq,
    /// `->`, `→`.
    Arrow,
    /// `(`.
    LParen,
    /// `)`.
    RParen,
    /// `[`.
    LBracket,
    /// `]`.
    RBracket,
    /// `{`.
    LBrace,
    /// `}`.
    RBrace,
    /// `,`.
    Comma,
    /// `;`.
    Semicolon,
    /// `√`.
    Sqrt,
    /// `∛`.
    Cbrt,
    /// `∜`.
    FourthRoot,

    /// A character that starts no token. Its error says why.
    Error,
    /// The end of the text, as the parser sees it; never in a tree.
    Eof,

    /// The whole text: its statements, separated by `;` or line breaks.
    Root,
    /// The statements of a block: a branch of `if`, the body of a loop or of
    /// `begin … end`.
    Body,
    /// An expression used as a statement.
    ExprStmt,
    /// `lhs := rhs`, maybe with a [`WhenClause`](SyntaxKind::WhenClause): an
    /// assignment or a definition, which lowering tells apart by the
    /// left-hand side (CAS-DOC-7 §4).
    Definition,
    /// `when guard`, after the right-hand side of a definition.
    WhenClause,
    /// `local a, b := e`.
    LocalStmt,
    /// One name of a `local` statement, with its optional `:= value`.
    LocalItem,
    /// `return`, with an optional value.
    ReturnStmt,
    /// `break`.
    BreakStmt,
    /// `continue`.
    ContinueStmt,
    /// `for name in e do … end`.
    ForStmt,
    /// `while c do … end`.
    WhileStmt,
    /// `params -> body`.
    Lambda,
    /// The parenthesised parameters of a lambda: `(x, y)`.
    ParamList,
    /// Operands with an operator between each two, all of one binding power
    /// (CAS-164): `or`; `and`; `+` and `-`; `*`, `/`, `.*` and `./`. So `a +
    /// b - c` is one node, and `a * b + c` a sum whose first operand is a
    /// product.
    BinaryExpr,
    /// `a = b`, an equation.
    Equation,
    /// `a < b <= c`: comparisons chained (`==`, `!=`, `<`, `<=`, `>`, `>=`,
    /// `in`).
    Comparison,
    /// `a..b`.
    RangeExpr,
    /// Operands side by side: `2x`, `x y`, `2(x + 1)`.
    ImplicitProduct,
    /// A prefix operator and its operand: `-x`, `√x`, `not p`.
    PrefixExpr,
    /// `b^e`, `b.^e`, `b²`.
    PowerExpr,
    /// `f(args)`: a call, or a product that lowering recognises
    /// (CAS-DOC-7 §3.1).
    Application,
    /// `v[i, j]`.
    IndexExpr,
    /// The delimited arguments of an application or an index.
    ArgList,
    /// `n!`.
    Factorial,
    /// `xs...`.
    Spread,
    /// `x :: type`.
    Typed,
    /// `(e)`.
    Paren,
    /// `[a, b]`.
    ListExpr,
    /// `{a, b}`.
    SetExpr,
    /// A number or a string.
    Literal,
    /// A name.
    Name,
    /// `if c then … elif c then … else … end`.
    IfExpr,
    /// `begin … end`.
    BlockExpr,
    /// Tokens the parser could not place, kept so that the tree stays
    /// lossless.
    ErrorNode,
}

impl SyntaxKind {
    /// Whitespace and comments. Line breaks are not trivia: the parser decides
    /// which of them separate statements.
    #[must_use]
    pub fn is_trivia(self) -> bool {
        matches!(
            self,
            SyntaxKind::Whitespace | SyntaxKind::LineComment | SyntaxKind::BlockComment
        )
    }

    /// The keywords, `if` to `minus`.
    #[must_use]
    pub fn is_keyword(self) -> bool {
        (SyntaxKind::IfKw..=SyntaxKind::MinusKw).contains(&self)
    }

    /// Whether this is the kind of a token rather than of a node.
    #[must_use]
    pub fn is_token(self) -> bool {
        self < SyntaxKind::Root
    }

    /// The keyword spelled `text`, if any. Keywords are ASCII and
    /// case-sensitive: `If` is a name.
    #[must_use]
    pub fn keyword(text: &str) -> Option<SyntaxKind> {
        Some(match text {
            "if" => SyntaxKind::IfKw,
            "then" => SyntaxKind::ThenKw,
            "elif" => SyntaxKind::ElifKw,
            "else" => SyntaxKind::ElseKw,
            "end" => SyntaxKind::EndKw,
            "for" => SyntaxKind::ForKw,
            "in" => SyntaxKind::InKw,
            "do" => SyntaxKind::DoKw,
            "while" => SyntaxKind::WhileKw,
            "begin" => SyntaxKind::BeginKw,
            "return" => SyntaxKind::ReturnKw,
            "break" => SyntaxKind::BreakKw,
            "continue" => SyntaxKind::ContinueKw,
            "local" => SyntaxKind::LocalKw,
            "when" => SyntaxKind::WhenKw,
            "and" => SyntaxKind::AndKw,
            "or" => SyntaxKind::OrKw,
            "not" => SyntaxKind::NotKw,
            "union" => SyntaxKind::UnionKw,
            "intersect" => SyntaxKind::IntersectKw,
            "minus" => SyntaxKind::MinusKw,
            _ => return None,
        })
    }

    /// How a token of this kind is written, for messages: `` `)` ``, `` `end` ``,
    /// or a description for tokens without a fixed text.
    #[must_use]
    pub fn describe(self) -> &'static str {
        match self {
            SyntaxKind::Number => "a number",
            SyntaxKind::String => "a string",
            SyntaxKind::Superscript => "a superscript",
            SyntaxKind::Ident => "a name",
            SyntaxKind::IfKw => "`if`",
            SyntaxKind::ThenKw => "`then`",
            SyntaxKind::ElifKw => "`elif`",
            SyntaxKind::ElseKw => "`else`",
            SyntaxKind::EndKw => "`end`",
            SyntaxKind::ForKw => "`for`",
            SyntaxKind::InKw => "`in`",
            SyntaxKind::DoKw => "`do`",
            SyntaxKind::WhileKw => "`while`",
            SyntaxKind::BeginKw => "`begin`",
            SyntaxKind::ReturnKw => "`return`",
            SyntaxKind::BreakKw => "`break`",
            SyntaxKind::ContinueKw => "`continue`",
            SyntaxKind::LocalKw => "`local`",
            SyntaxKind::WhenKw => "`when`",
            SyntaxKind::AndKw => "`and`",
            SyntaxKind::OrKw => "`or`",
            SyntaxKind::NotKw => "`not`",
            SyntaxKind::UnionKw => "`union`",
            SyntaxKind::IntersectKw => "`intersect`",
            SyntaxKind::MinusKw => "`minus`",
            SyntaxKind::Plus => "`+`",
            SyntaxKind::Minus => "`-`",
            SyntaxKind::Star => "`*`",
            SyntaxKind::Slash => "`/`",
            SyntaxKind::Caret => "`^`",
            SyntaxKind::DotStar => "`.*`",
            SyntaxKind::DotSlash => "`./`",
            SyntaxKind::DotCaret => "`.^`",
            SyntaxKind::Bang => "`!`",
            SyntaxKind::Eq => "`=`",
            SyntaxKind::EqEq => "`==`",
            SyntaxKind::Neq => "`!=`",
            SyntaxKind::Lt => "`<`",
            SyntaxKind::Le => "`<=`",
            SyntaxKind::Gt => "`>`",
            SyntaxKind::Ge => "`>=`",
            SyntaxKind::DotDot => "`..`",
            SyntaxKind::Ellipsis => "`...`",
            SyntaxKind::ColonColon => "`::`",
            SyntaxKind::ColonEq => "`:=`",
            SyntaxKind::Arrow => "`->`",
            SyntaxKind::LParen => "`(`",
            SyntaxKind::RParen => "`)`",
            SyntaxKind::LBracket => "`[`",
            SyntaxKind::RBracket => "`]`",
            SyntaxKind::LBrace => "`{`",
            SyntaxKind::RBrace => "`}`",
            SyntaxKind::Comma => "`,`",
            SyntaxKind::Semicolon => "`;`",
            SyntaxKind::Sqrt => "`√`",
            SyntaxKind::Cbrt => "`∛`",
            SyntaxKind::FourthRoot => "`∜`",
            SyntaxKind::Eof => "the end of the text",
            SyntaxKind::Newline => "a line break",
            SyntaxKind::Whitespace
            | SyntaxKind::LineComment
            | SyntaxKind::BlockComment
            | SyntaxKind::Error => "a token",
            _ => "a node",
        }
    }
}

/// The language of the syntax trees, for `rowan`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum CasLanguage {}

impl rowan::Language for CasLanguage {
    type Kind = SyntaxKind;

    fn kind_from_raw(raw: rowan::SyntaxKind) -> SyntaxKind {
        SyntaxKind::ALL[usize::from(raw.0)]
    }

    fn kind_to_raw(kind: SyntaxKind) -> rowan::SyntaxKind {
        rowan::SyntaxKind(kind as u16)
    }
}

/// A node of a syntax tree, with its position in the text.
pub type SyntaxNode = rowan::SyntaxNode<CasLanguage>;
/// A token of a syntax tree, with its position in the text.
pub type SyntaxToken = rowan::SyntaxToken<CasLanguage>;
/// A node or a token.
pub type SyntaxElement = rowan::SyntaxElement<CasLanguage>;

impl From<SyntaxKind> for rowan::SyntaxKind {
    fn from(kind: SyntaxKind) -> rowan::SyntaxKind {
        <CasLanguage as rowan::Language>::kind_to_raw(kind)
    }
}
