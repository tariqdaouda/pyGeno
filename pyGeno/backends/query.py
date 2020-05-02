class Filter:
    """docstring for Filter"""
    def __init__(self, statement, operator=None, right_statement=None):
        super(Filter, self).__init__()
        self.statement = statement
        self.operator = operator
        self.right_statement = right_statement

    def _mix(self, statement, operator):
        if not isinstance(statement, Filter):
            raise ValueError("statement must be of type Filter, got: '%s'" % statement)
        return Filter(self, operator, statement)

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

f = Filter

if __name__ == '__main__':
   q =  ( f("a > 1") & f("b < 6") ) | ( f("d = 1") & f("c != 6") ) & f("j = 0") 
   print(q)
