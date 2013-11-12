#!/usr/bin/env python

# test the fitting methods.
from snpy import get_sn,model
import types

models = []
for item in model.__dict__:
   obj = model.__dict__[item]
   if type(obj) is types.ClassType:
      if issubclass(obj, model.model):
         models.append(item)

for model in models:
   if model == 'model': continue
   print "*******************",model,"******************"
   s = get_sn('SN2006ax_kcorr.snpy')
   if model == 'EBV_model': 
      stype='dm15'
   else:
      stype='st'
   s.choose_model(model, stype=stype)
   s.fit(['u','B','V','g','r','i','Y','J','H'], kcorr=False)
   s.summary()
   s.save('SN2006ax_%s.snpy' % model)
