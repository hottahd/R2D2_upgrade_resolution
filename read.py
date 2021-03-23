import numpy as np
import f90nml

param = f90nml.read('namelist.in')
caseid = param['info']['caseid_out']
margin = 2

ixg = param['upgrade']['ix00'] + 2*margin
jxg = param['upgrade']['jx00'] + 2*margin
kxg = param['upgrade']['kx00'] + 2*margin
mtype = 9

f = open('../run/'+caseid+'/data/qq/qq.dac.e','rb')
qq = np.fromfile(f,'<d',mtype*ixg*jxg*kxg).reshape((ixg,jxg,kxg,mtype),order='F')
