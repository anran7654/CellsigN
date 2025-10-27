import re

s = "a,b\tc,d\te"
parts = re.split('[\t,]', s)
print(parts)
