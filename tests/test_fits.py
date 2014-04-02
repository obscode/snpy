#!/usr/bin/env python

# test the fitting methods.
import matplotlib
matplotlib.use('Agg')
from snpy import get_sn,model
import types

models = []
for item in model.__dict__:
   obj = model.__dict__[item]
   if type(obj) is types.ClassType:
      if issubclass(obj, model.model):
         models.append(item)

for model in models:
   if model in ['model','Rv_model','color_model']: continue
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
   
   if 'Bs' in s.model.rbs:
      s.restbands['B'] = 'Bs'
      s.restbands['V'] = 'Vs'
      s.restbands['r'] = 'Rs'
      s.restbands['i'] = 'Is'
   s.fit(['u','B','V','g','r','i','Y','J','H'], kcorr=False)
   s.summary()
   s.save('SN2006ax_%s_st.snpy' % model)

