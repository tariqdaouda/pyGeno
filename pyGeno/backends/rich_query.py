class RichFilter:
    """docstring for RichFilter"""
    def __init__(self, statement, operator=None, right_statement=None):
        super(RichFilter, self).__init__()
        self.statement = statement
        self.operator = operator
        self.right_statement = right_statement

    def _mix(self, statement, operator):
        if not isinstance(statement, RichFilter):
            raise ValueError("statement must be of type RichFilter, got: '%s'" % statement)
        return RichFilter(self, operator, statement)

    def __and__(self, statement):
        return self._mix(statement, "and")

    def __or__(self, statement):
        return self._mix(statement, "or")
    
    def to_str(self):
        if self.operator:
            return "(%s %s %s)" % (str(self.statement), str(self.operator), str(self.right_statement)) 
        return str(self.statement)

    def __repr__(self):
        return self.to_str()

_f = RichFilter

if __name__ == '__main__':
   q =  ( _f("a > 1") & _f("b < 6") ) | ( _f("d = 1") & _f("c != 6") ) & _f("j = 0") 
   print(q)
