#functions to automatcly generate SQL statements

def search(dataType, dictionary = None, suffix = '') :
	if dictionary != None :
		conditionsStr = ''
		conditionsValues = []
		for k in dictionary :
			value = dictionary[k].strip()
			if dictionary[k].find('...') > -1 :
				operator = "LIKE"
				value = value.replace('...', '%')
			elif value[0] == '>' or value[0] == '<' or value[0] == '=' :
				operator = value[0]
				value = value[1:]
			elif value[2:] == '>=' or value[2:] == '<=' :
				operator = value[2:]
				value = value[2:]
			else :
				operator = '='
			if value[:3].lower() == 'or ' or value[:4].lower() == 'and ' :
				conditionsStr = '%s %s ? ' % (k, operator)
				conditionsValues.append(value)
			else :
				conditionsStr += '%s %s ? ' % (k, operator)
				conditionsValues.append(value)
		
		sql = 'SELECT * FROM %s WHERE %s %s;' % (dataType, conditionsStr, suffix)
		return (sql, conditionsValues)
	else :
		return ('SELECT * FROM %s %s;' % suffix, ())
	
if __name__ == '__main__' :
	print search('gene', {'symbol' : 'IG...', 'and length' : '< 200', 'or chromosome' : '3'})
