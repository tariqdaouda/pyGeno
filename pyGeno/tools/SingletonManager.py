#This thing is wonderful

objects = {}
def add(obj, objName='') :
	
	if objName == '' :
		key = obj.name
	else :
		key = objName
		
	if key not in objects :
		objects[key] = obj

	return obj

def contains(k) :
	return k in objects
	
def get(objName) :
	try :
		return objects[objName]
	except :
		return None
