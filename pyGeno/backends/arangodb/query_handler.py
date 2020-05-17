from ..backend_abs import QueryHandler_ABS

class QueryHandler(QueryHandler_ABS):
    """
    Handles different query types
    """
    def __init__(self):
        pass

    def rich_query(self, pygeno_filter, object_name):
        def _rec_rename(pygeno_filter, object_name):
            if pygeno_filter.right_statement:
                _rec_rename(pygeno_filter.statement, object_name)
                _rec_rename(pygeno_filter.right_statement, object_name)
            else :
                pygeno_filter.statement = "%s.%s" %(object_name, pygeno_filter.statement)

        _rec_rename(pygeno_filter, object_name)
        fi = pygeno_filter.to_str().replace("and", "&&").replace("or", "||")
        aql_filters = "FILTER %s" % fi

        aql = """FOR {obj} IN {collection}
        {filters}
        RETURN {obj} 
        """.format(collection=self.collection_name, obj=object_name, filters=aql_filters)
    
        objects = self.db.AQLQuery(aql, batchSize=1000, rawResults=True)
        for obj in objects;
            yield obj_type(obj)

    def lazy_query(self, anchor_type, type_name, **lazy_args):
        return dict_query(anchor_type, type_name, lazy_args)

    def dict_query(self, anchor_type, type_name, dct):
        def _find_operator(key):
            for op in [">", ">=", "<". "<=", "=="]:
                if op in key :
                    return op
            return "=="

        filters = []
        for key, value in dct.items():
            op = _find_operator(key)
            filters.append(
                "FILTER obj2.{field} {op} {value}".format(format=key, op=op, value=value)
            )
        aql_filters = "\n".join(filters)
        aql == """
        FOR obj1 IN {anchor_col}
            FOR link IN {anchor_col}_to_{other_col}
                FILTER link._from == obj1._id OR link._to == obj1._id       
                FOR obj2 IN {other_col}
                    FILTER link._from == obj2._id OR link._to == obj2._id
                    {filters}
                    RETURN obj2 
        """.format(anchor_col=anchor_type, other_col=type_name, filters=aql_filters)

        objects = self.db.AQLQuery(aql, batchSize=1000, rawResults=True)
        for obj in objects;
            yield obj_type(obj)
