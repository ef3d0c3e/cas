//! Token kinds. The parser (CAS-118) adds the kinds of tree nodes.

/// The kind of a token.
///
/// Aliases share a kind with what they stand for: `×` is a [`Star`] and `≤` a
/// [`Le`], with their own text (CAS-DOC-7 §1). `π`, `ℯ` and `∞` are
/// [`Ident`]s of one character; [`name`](crate::name) gives `pi`, `e` and
/// `inf` for them.
///
/// [`Star`]: SyntaxKind::Star
/// [`Le`]: SyntaxKind::Le
/// [`Ident`]: SyntaxKind::Ident
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u16)]
pub enum SyntaxKind {
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

    /// The keywords, `if` to `not`.
    #[must_use]
    pub fn is_keyword(self) -> bool {
        (SyntaxKind::IfKw..=SyntaxKind::NotKw).contains(&self)
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
            _ => return None,
        })
    }
}
