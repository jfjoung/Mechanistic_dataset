__all__ = ['ReactionError', 'NoReactantError', 'NoAcidBaseError', 'TooManyTemplates']

class ReactionError(Exception):
    """Base class for exceptions in this module."""
    pass

class NoReactantError(ReactionError):
    """Exception raised when there are no reactants in the flask."""
    pass

class NoAcidBaseError(ReactionError):
    """Exception raised when there are no acids nor bases in the flask."""
    pass

class TooManyTemplates(ReactionError):
    """Exception raised when there are too many possible templates after enumeration."""
    pass
