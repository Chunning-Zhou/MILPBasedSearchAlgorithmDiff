from ciphers.present import class_present
from ciphers.gift import class_gift
from ciphers.rectangle import class_rectangle
from ciphers.lblock import class_lblock
from ciphers.twine import class_twine
''''
param cipher: a bloch cipher class
'''

# cipher =  class_present()
# cipher =  class_gift()
# cipher =  class_rectangle()
cipher = class_lblock()
# cipher = class_twine()

'''
param goal:  
goal = AS: calculate the minimum number of act
goal = DC: search for the best differential characteriatic
'''
goal = "DC"


