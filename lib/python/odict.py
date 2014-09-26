#!/task/imd/python/bin/python
#
# odict.py -- Version=2.1
#             July/22/2008 YoungMu Jeon
#             mailto: jeon@fusion.gat.com
#
"""
A generalized, standard "ordered" dictionary

OPTIONS:
 - CASE   : if None   , case-sensitive (default)
            if "upper", case-insensitive and upper-character displayed
            if "lower", case-insensitive and lower-character displayed
"""
"""
UPDATE HISTORY:
- A method 'copy()' upgraded in order to copy real values, not references.
                                                                  - 20080722 YMJ
- New option for CASE-SENSITIVITY added and some basic corrections made
                                                                  - 20080412 YMJ
- First created                                                   - 20080201 YMJ
"""

class odict:
   """ """
   def __init__(self,dict=None,case=None,**kwargs):
       self.data = {}
       self.okeys= []
       self.case = case
       if dict:
          self.update(dict)
       if len(kwargs):
          self.update(kwargs)
   
   def __repr__(self):
       str = "{"
       for key in self.okeys:
           str += "'%s': %s, " %(key,self[key])
       if  len(str)==1: str += "}"
       else: str = str[:-2]+"}"
       
       return str
   
   def __len__(self): return len(self.data)

   def __getitem__(self, key):
       if self.case:
          if self.case.lower()=="upper": key = key.upper()
          if self.case.lower()=="lower": key = key.lower()
       if key in self.okeys: return self.data[key]
       else: print "key not found"; return

   def __setitem__(self, key, item):
       if  self.case:
           if self.case.lower()=="upper": key = key.upper()
           if self.case.lower()=="lower": key = key.lower()
       if  hasattr(item,"items"):
           self.data[key] = odict(item,case=self.case)
       else:
           self.data[key] = item
       
       if not key in self.okeys: self.okeys.append(key)
   
   def __delitem__(self, key):
       if self.case:
          if self.case.lower()=="upper": key = key.upper()
          if self.case.lower()=="lower": key = key.lower()
       if key in self.okeys:
          self.okeys.remove(key)
          del self.data[key]
   
   def __add__(self,other):
       new_odict = self.copy()
       new_odict.update(other)
       return new_odict
   
   def clear(self): self.data.clear(); self.okeys=[]
   def copy(self):
       """
       This copy is not same with intrinsic copy() of python. Intrinsic copy() just copies
       its value, so that for a dictionary or a list it copies its reference, not value. 
       Here this copy() is designed to copy its corresponding value. 
       Note that in most cases it works well, but if it has a dictionary as its value, then
       it copies its reference. (It's impossible to generalize this functionality to a
       dictionary-type value, because there is no way to figure out how many levels that
       dictionary has.)
       """
       import copy
       return copy.deepcopy(self)

   def keys(self): return self.okeys[:]

   def items(self):
       item = []
       for k in self.okeys: item.append((k,self[k]))
       return item

   def has_key(self, key):
       if self.case:
          if self.case.lower()=="upper": key = key.upper()
          if self.case.lower()=="lower": key = key.lower()
       if key in self.okeys: return True
       return False

   def update(self, dict=None, **kwargs):
       """ This update means adding new ones and replacing existing ones. """
       # For dict
       if dict:
          for k,v in dict.items():
             if self.case:
                if self.case.lower()=="upper": k = k.upper()
                if self.case.lower()=="lower": k = k.lower()
             if hasattr(v,"items"): self.data[k] = odict(v,case=self.case)
             else                 : self.data[k] = v
             if not k in self.okeys: self.okeys.append(k)
       # For kwargs
       if len(kwargs):
          for k, v in kwargs.items():
             if self.case:
                if self.case.lower()=='upper': k = k.upper()
                if self.case.lower()=='lower': k = k.lower()
             if hasattr(v,"items"): self.data[k] = odict(v,case=self.case)
             else                 : self.data[k] = v
             if not k in self.okeys: self.okeys.append(k)

   def get(self,key,notkey=None):
       """
       Get the value corresponding to the key. If the given key is not
       recognized, then the notkey's value is returned
       """
       if self.case:
          if self.case.lower()=="upper": key = key.upper()
          if self.case.lower()=="lower": key = key.lower()
       if   key in self.okeys: return self[key]
       elif notkey: return notkey
       else: return None

   def pop(self,key,notkey=None):
       if self.case:
          if self.case.lower()=="upper": key = key.upper()
          if self.case.lower()=="lower": key = key.lower()
       if key in self.okeys:
          self.okeys.remove(key)
          return self.data.pop(key)
       elif notkey: return notkey
       else: return None
   
   def popitem(self):
       if len(self.okeys):
          key = self.okeys.pop()
          val = self.data.pop(key)
          return (key,val)
       else: return (None,None)
   
   def __contains__(self, key):
       if self.case:
          if self.case.lower()=="upper": key = key.upper()
          if self.case.lower()=="lower": key = key.lower()
       return self.data._contains_(key)

