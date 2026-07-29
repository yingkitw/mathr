use thiserror::Error;

pub type Result<T> = std::result::Result<T, MathError>;

#[derive(Debug, Error)]
pub enum MathError {
    #[error("parse error: {0}")]
    Parse(String),

    #[error("evaluation error: {0}")]
    Eval(String),

    #[error("unknown variable: {0}")]
    UnknownVariable(String),

    #[error("unknown function: {0}")]
    UnknownFunction(String),

    #[error("domain error: {0}")]
    Domain(String),

    #[error("non-convergent: {0}")]
    NotConvergent(String),

    #[error("invalid argument: {0}")]
    InvalidArgument(String),

    #[error("io error: {0}")]
    Io(#[from] std::io::Error),

    #[error("plot error: {0}")]
    Plot(String),

    #[error("{0}")]
    Other(String),
}

impl From<anyhow::Error> for MathError {
    fn from(e: anyhow::Error) -> Self {
        MathError::Other(e.to_string())
    }
}