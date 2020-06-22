from ..backend_abs import QueryHandler_ABS

class QueryHandler(QueryHandler_ABS):
    """
    Handles different query types
    """

    def __init__(self, *args, **kwargs):
        super(QueryHandler, self).__init__(*args, **kwargs)
        self.db = self.database_configuration.database

    def get_example(self, object_type):
        return self.db[object_type].fetchFirstExample({}).result

    def lazy_query(self, anchor_type, fetch_type=None, **lazy_args):
        return self.dict_query(anchor_type, fetch_type, lazy_args)

    def _find_operator(sels, key):
        for op in [">=", "<=", "==", ">", "<"]:
            if op in key :
                return op
        return "=="

    def _dct_to_aql_filters(self, dct=None, obj_name="obj"):
        if dct is None :
            return ""

        filters = []
        for key, value in dct.items():
            op = self._find_operator(key)
            if type(value) is str :
                value = '"%s"' % value
            key = key.replace(op, "")
            filters.append(
                "FILTER {obj}.{field} {op} {value}".format(field=key, op=op, value=value, obj=obj_name)
            )
        aql_filters = "\n".join(filters)
        return aql_filters

    def _get_dict_result(self, anchor_type, dct=None):
        aql_filters = self._dct_to_aql_filters(dct, obj_name="obj")

        aql = """
        FOR obj IN {anchor_col}
            {filters}
            RETURN obj 
        """.format(anchor_col=anchor_type, filters=aql_filters)

        objects = self.db.AQLQuery(aql, batchSize=1000, rawResults=True)
        for obj in objects:
            yield obj

    def _get_join_result(self, anchor_type, fetch_type, aql_filters):
        link_col_name = self.database_configuration.get_link_collection_name(anchor_type, fetch_type)
        aql = """
        FOR obj1 IN {anchor_col}
            FOR link IN {link_col_name}
                FILTER link._from == obj1._id OR link._to == obj1._id       
                FOR obj2 IN {other_col}
                    FILTER link._from == obj2._id OR link._to == obj2._id
                    {filters}
                    RETURN obj2 
        """.format(anchor_col=anchor_type, link_col_name=link_col_name,other_col=fetch_type, filters=aql_filters)

        objects = self.db.AQLQuery(aql, batchSize=1000, rawResults=True)
        for obj in objects:
            yield obj
    
    def _get_dict_join_result(self, anchor_type, fetch_type, dct=None):
        aql_filters = self._dct_to_aql_filters(dct, obj_name="obj2")
        return self._get_join_result(anchor_type, fetch_type, aql_filters)

    def dict_query(self, anchor_type, fetch_type=None, dct=None,):
        if fetch_type is None :
            return self._get_dict_result(anchor_type, dct)
        else:
            return self._get_dict_join_result(anchor_type, fetch_type, dct)

    def _get_rich_query_filters(self, pygeno_filter, obj_name):
        def _rec_rename(pygeno_filter, obj_name):
            if pygeno_filter.right_statement:
                _rec_rename(pygeno_filter.statement, obj_name)
                _rec_rename(pygeno_filter.right_statement, obj_name)
            else :
                pygeno_filter.statement = "%s.%s" %(obj_name, pygeno_filter.statement)

        _rec_rename(pygeno_filter, obj_name)
        fi = pygeno_filter.to_str().replace("and", "&&").replace("or", "||")
        aql_filters = "FILTER %s" % fi
        return aql_filters

    def rich_query(self, anchor_type, fetch_type, pygeno_filter):
        aql_filters = self._get_rich_query_filters(pygeno_filter, obj_name='obj2')
        return self._get_join_result(anchor_type, fetch_type, aql_filters)
