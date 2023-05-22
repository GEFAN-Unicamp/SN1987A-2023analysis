import numpy as np
from flavor_conversion import *
from scipy import interpolate


#Function that do a marginalized profile over a parameter "param" given a Minuit object "m"
def param_profile(m,param,min_val,max_val,size=50,iterate=5,subtract_min=True):
    values=np.linspace(min_val,max_val,size)
    like=[]
    values_aux=[]
    other_params=[]
    
    intial_values=np.array(m.values)
    m.fixed[param] = True
    
    below=values[values<=m.values[param]]
    above=values[values>m.values[param]]
    
    #Scaning below de BF value
    for value in below[::-1]:
        m.values[param]=value
        m.migrad(iterate = iterate)
        if m.valid:
            like.append(m.fval)
            values_aux.append(value)
            other_params.append(np.array(m.values))
        else:
            print("Faield to converge in ",param," = ",value)

    
    #Reversing the arrays
    like=like[::-1]
    values_aux=values_aux[::-1]
    other_params=other_params[::-1]
    #Getting back to the initial BF values
    m.values=intial_values
    
    #Scaning above de BF value
    for value in above:
        m.values[param]=value
        m.migrad(iterate = iterate)
        if m.valid:
            like.append(m.fval)
            values_aux.append(value)
            other_params.append(np.array(m.values))
        else:
            print("Faield to converge in ",param," = ",value)
            
    m.fixed[param] = False
    m.values=intial_values
    
    like=np.array(like)
    if subtract_min==True:
        like=like-np.amin(like)
    
    return values_aux,like,other_params

def param_profile_mixed(m,m1,m2,param,min_val,max_val,size=50,iterate=5,subtract_min=True):
    values=np.linspace(min_val,max_val,size)
    like=[]
    values_aux=[]
    other_params=[]
    
    intial_values=np.array(m.values)
    intial_values1=np.array(m1.values)
    intial_values2=np.array(m2.values)
    m.fixed[param] = True
    m1.fixed[param] = True
    m2.fixed[param] = True
    
    below=values[values<=m.values[param]]
    above=values[values>m.values[param]]
    
    #Scaning below de BF value
    for value in below[::-1]:
        m1.values[param]=value
        m1.migrad(iterate = iterate)
        m2.values[param]=value
        m2.migrad(iterate=iterate)
        
        if m1.fval<m2.fval:
            aux=list(m1.values)
            aux.append(U[0][0]**2)
            m.values=aux
        else:
            aux=list(m2.values)
            aux.append(U[0][2]**2)
            m.values=aux
            
        m.values[param]=value
        m.migrad(iterate = iterate)    
        
        if m.valid:
            like.append(m.fval)
            values_aux.append(value)
            other_params.append(np.array(m.values))
        else:
#             like.append(m.fval)
#             values_aux.append(value)
#             other_params.append(np.array(m.values))
            print("Faield to converge in ",param," = ",value)

    
    #Reversing the arrays
    like=like[::-1]
    values_aux=values_aux[::-1]
    other_params=other_params[::-1]
    #Getting back to the initial BF values
    m.values=intial_values
    m1.values=intial_values1
    m2.values=intial_values2
    
    #Scaning above de BF value
    for value in above:
        m1.values[param]=value
        m1.migrad(iterate = iterate)
        m2.values[param]=value
        m2.migrad(iterate=iterate)
        
        if m1.fval<m2.fval:
            aux=list(m1.values)
            aux.append(U[0][0]**2)
            m.values=aux
        else:
            aux=list(m2.values)
            aux.append(U[0][2]**2)
            m.values=aux
            
        m.values[param]=value
        m.migrad(iterate = iterate) 
        
        if m.valid:
            like.append(m.fval)
            values_aux.append(value)
            other_params.append(np.array(m.values))
        else:
#             like.append(m.fval)
#             values_aux.append(value)
#             other_params.append(np.array(m.values))
            print("Faield to converge in ",param," = ",value)
            
    m.fixed[param] = False
    m1.fixed[param] = False
    m2.fixed[param] = False
    m.values=intial_values
    m1.values=intial_values1
    m2.values=intial_values2
    
    like=np.array(like)
    if subtract_min==True:
        like=like-np.amin(like)
    
    return values_aux,like,other_params

def param_profile_v3(m,param,min_val,max_val,size=50,iterate=5,subtract_min=True):
    values=np.linspace(min_val,max_val,size)
    like=[]
    values_aux=[]
    other_params=[]
    
    intial_values=np.array(m.values)
    m.fixed[param] = True
    
    below=values[values<=m.values[param]]
    above=values[values>m.values[param]]
    
    #Scaning below de BF value
    for value in below[::-1]:
        m.values[param]=value
        #t_off_I=0.5s
        m.fixed["x8"] = True
        m.values["x8"]=0.5
        m.migrad(iterate = iterate//2)
        values_1=list(m.values)
        fval_1=m.fval
        
        #t_off_I=0.0s
        m.fixed["x8"] = True
        m.values["x8"]=0
        m.migrad(iterate = iterate//2)
        values_2=list(m.values)
        fval_2=m.fval
        
        #Best
        m.fixed["x8"] = False
        if fval_1<=fval_2:
            m.values=values_1
            m.migrad(iterate = iterate)
        else:
            m.values=values_2
            m.migrad(iterate = iterate)
        
        like.append(m.fval)
        values_aux.append(value)
        other_params.append(np.array(m.values))
        
        if not m.valid:
            print("Faield to converge in ",param," = ",value)

    
    #Reversing the arrays
    like=like[::-1]
    values_aux=values_aux[::-1]
    other_params=other_params[::-1]
    #Getting back to the initial BF values
    m.values=intial_values
    
    #Scaning above de BF value
    for value in above:
        m.values[param]=value
        #t_off_I=0.5s
        m.fixed["x8"] = True
        m.values["x8"]=0.5
        m.migrad(iterate = iterate//2)
        values_1=list(m.values)
        fval_1=m.fval
        
        #t_off_I=0.0s
        m.fixed["x8"] = True
        m.values["x8"]=0
        m.migrad(iterate = iterate//2)
        values_2=list(m.values)
        fval_2=m.fval
        
        #Best
        m.fixed["x8"] = False
        m.values[param]=value
        if fval_1<=fval_2:
            m.values=values_1
            m.migrad(iterate = iterate)
        else:
            m.values=values_2
            m.migrad(iterate = iterate)
            
    m.fixed[param] = False
    m.values=intial_values
    
    like=np.array(like)
    if subtract_min==True:
        like=like-np.amin(like)
    
    return values_aux,like,other_params



################ Algorithm to make a 2D path starting from the BF foint ###########
def ordering(ni,nj,iini,jini):
    
  ltmp=np.zeros((ni+2,nj+2),dtype=int)

  #Defining borders?
  for i in range(ni+2):
    ltmp[i][0]=1
    ltmp[i][nj+1]=1
  for j in range(nj+2):
    ltmp[0][j]=1
    ltmp[ni+1][j]=1
  
  k=[]
  i=iini+1
  j=jini+1
  k.append([i-1,j-1])
  ltmp[i][j]=1  #Is it 1 for already explored point?
  
  # step[new direction, old direction], counterclock
  step1=[[[0,1],[1,0]],[[-1,0],[0,1]],[[0,-1],[-1,0]],[[1,0],[0,-1]]]
  # step[new direction, old direction], clock
  step2=[[[0,-1],[1,0]],[[-1,0],[0,-1]],[[0,1],[-1,0]],[[1,0],[0,1]]]

  if iini<ni//2 and nj//2<jini:
#     print("left up")
    step=step2.copy()
    id=3
  elif iini>ni//2 and nj//2<jini:
#     print("right up")
    step=step1.copy()
    id=1
  elif iini<ni//2 and jini<nj//2:
    step=step1.copy()
    id=3
  else:
    id=0
    step=step1.copy()
    
    #   print(ni//2,iini,nj//2,jini)
  if iini<ni//2 and nj//2<jini:
#     print("left up")
    step=step2.copy()
    id=3
  elif iini>ni//2 and nj//2<jini:
#     print("right up")
    step=step1.copy()
    id=1
  elif iini<ni//2 and jini<nj//2:
    step=step1.copy()
    id=3
  else:
    id=0
    step=step1.copy()

  for ntmp in range(0,ni*nj):
    itst=0
    # changing direction
    ii,jj=step[id][0][0],step[id][0][1]
    if(ltmp[i+ii][j+jj]==0):
      id=id+1
      if(id==4):id=0
      # keeping direction
    else:
      ii,jj=step[id][1][0],step[id][1][1]
      # inverting direction
      if(ltmp[i+ii][j+jj]==1):
        if(step == step1): 
          step=step2.copy()
        elif(step == step2): 
          step=step1.copy()
        itst=1
        for itmp in range(0,4):
          ii,jj=step[itmp][0][0],step[itmp][0][1]
          if(ltmp[i+ii][j+jj]==0):
            id=itmp
            itst=0
        ii,jj=step[id][0][0],step[id][0][1]
        id=id+1
        if(id==4):id=0
            
    if(itst==1):
      if(ntmp != ni*nj-1):
#         print
#         print('breaking at:', ntmp)
        for kk in range(ni):
          for kl in range(nj):
            if(ltmp[kk+1][kl+1] == 0):k.append([kk,kl])
#       else:
#         print
#         print('complete!')
      break    

    i=i+ii
    j=j+jj        
    ltmp[i][j]=1
    k.append([i-1,j-1])
    ltmp[i][j]=1
  
  return k


#Function that do a marginalized 2D scan over the parameters "param1" and "param2" given a Minuit object "m"
def two_dim_scan(m,info1,info2,iterate=5,subtract_min=True,n_interp=100):
    param1,min_val1,max_val1,size1=info1
    param2,min_val2,max_val2,size2=info2
    
    intial_values=np.array(m.values)
    x_BF,y_BF=intial_values[param1],intial_values[param2]
    
    
    #Method 2  
    values1=np.linspace(min_val1,max_val1,size1)
    values2=np.linspace(min_val2,max_val2,size2)
    x_i=(np.abs(values1 - x_BF)).argmin()
    y_i=(np.abs(values2 - y_BF)).argmin()
    k=ordering(size1,size2,x_i,y_i)
    
    
    m.fixed[param1] = True
    m.fixed[param2] = True
    like=[]
    
    for i,j in k:
        m.values[param1]=values1[i]
        m.values[param2]=values2[j]
        m.migrad(iterate=iterate)
        like.append([values1[i],values2[j],m.fval,m.valid])
        if not m.valid:
            print("Faield to converge in ",param1,",",param2," = ",values1[i],",",values2[j])
    m.fixed[param1] = False
    m.fixed[param2] = False
    m.values=intial_values
    
    L=np.array(like)
    
    #Getting unique values
    x_aux=np.unique(L[:,0])
    y_aux=np.unique(L[:,1])
    L_aux=[]
    flag_aux=[]
    #Ordering
    for x_i in x_aux:
        L_aux2=[]
        flag_aux2=[]
        for y_i in y_aux:
            L_aux2.append(L[np.logical_and(L[:,0]==x_i,L[:,1]==y_i),2][0])
            flag_aux2.append(L[np.logical_and(L[:,0]==x_i,L[:,1]==y_i),3][0])
        L_aux.append(L_aux2)
        flag_aux.append(flag_aux2)
    
    #Interpolation

    f = interpolate.interp2d(x_aux, y_aux, L_aux, kind='cubic')    
    x_new = np.linspace(x_aux[0], x_aux[-1], n_interp)
    y_new = np.linspace(y_aux[0], y_aux[-1], n_interp)
    L_interp = f(x_new, y_new)
    
    if subtract_min==True:
        L_interp=L_interp-np.min(L_interp)
        L_aux=L_aux-np.min(L_aux)
    
    return [x_aux,x_new],[y_aux,y_new],[L_aux,L_interp], k, flag_aux


def two_dim_scan_mixed(m,m1,m2,info1,info2,iterate=5,subtract_min=True,n_interp=100):
    param1,min_val1,max_val1,size1=info1
    param2,min_val2,max_val2,size2=info2
    
    intial_values1=np.array(m1.values)
    intial_values2=np.array(m2.values)
    intial_values=np.array(m.values)
    x_BF,y_BF=intial_values[param1],intial_values[param2]
    
    
    #Method 2  
    values1=np.linspace(min_val1,max_val1,size1)
    values2=np.linspace(min_val2,max_val2,size2)
    x_i=(np.abs(values1 - x_BF)).argmin()
    y_i=(np.abs(values2 - y_BF)).argmin()
    k=ordering(size1,size2,x_i,y_i)
    
    
    m1.fixed[param1] = True
    m1.fixed[param2] = True
    m2.fixed[param1] = True
    m2.fixed[param2] = True
    m.fixed[param1] = True
    m.fixed[param2] = True
    like=[]
    
    for i,j in k:
        m1.values[param1]=values1[i]
        m1.values[param2]=values2[j]
        m1.migrad(iterate=iterate)
        m2.values[param1]=values1[i]
        m2.values[param2]=values2[j]
        m2.migrad(iterate=iterate)
        
        if m1.fval<m2.fval:
            aux=list(m1.values)
            aux.append(U[0][0]**2)
            m.values=aux
        else:
            aux=list(m2.values)
            aux.append(U[0][2]**2)
            m.values=aux
            
        m.migrad(iterate=iterate)
        like.append([values1[i],values2[j],m.fval,m.valid])
        if not m.valid:
            print("Faield to converge in ",param1,",",param2," = ",values1[i],",",values2[j])
            
    m1.fixed[param1] = False
    m1.fixed[param2] = False
    m1.values=intial_values1
    m2.fixed[param1] = False
    m2.fixed[param2] = False
    m2.values=intial_values2
    m.fixed[param1] = False
    m.fixed[param2] = False
    m.values=intial_values
    
    L=np.array(like)
    
    #Getting unique values
    x_aux=np.unique(L[:,0])
    y_aux=np.unique(L[:,1])
    L_aux=[]
    flag_aux=[]
    #Ordering
    for x_i in x_aux:
        L_aux2=[]
        flag_aux2=[]
        for y_i in y_aux:
            L_aux2.append(L[np.logical_and(L[:,0]==x_i,L[:,1]==y_i),2][0])
            flag_aux2.append(L[np.logical_and(L[:,0]==x_i,L[:,1]==y_i),3][0])
        L_aux.append(L_aux2)
        flag_aux.append(flag_aux2)
    
    #Interpolation

    f = interpolate.interp2d(x_aux, y_aux, L_aux, kind='cubic')    
    x_new = np.linspace(x_aux[0], x_aux[-1], n_interp)
    y_new = np.linspace(y_aux[0], y_aux[-1], n_interp)
    L_interp = f(x_new, y_new)
    
    if subtract_min==True:
        L_interp=L_interp-np.min(L_interp)
        L_aux=L_aux-np.min(L_aux)
    
    return [x_aux,x_new],[y_aux,y_new],[L_aux,L_interp], k, flag_aux


####################################### Old Functions ##############################
# ## Function to do a scan over 2 parameters minimizing the others
# def two_dim_scan(m,info1,info2):
#     param1,min_val1,max_val1,size1=info1
#     param2,min_val2,max_val2,size2=info2
    
#     values1=np.linspace(min_val1,max_val1,size1)
#     values2=np.linspace(min_val2,max_val2,size2)
    
#     intial_values=np.array(m.values)
#     like=[]

#     m.fixed[param1] = True
#     m.fixed[param2] = True
    
#     for i in range(len(values1)):
#         like.append([])
#         for j in range(len(values2)):
#             m.values=intial_values
#             m.values[param1]=values1[i]
#             m.values[param2]=values2[j]
#             m.migrad()
#             if not m.valid:
#                 print("Faield to converge in ",param1,",",param2," = ",values1[i],",",values2[j])
#             like[i].append(m.fval)
#     m.fixed[param1] = False
#     m.fixed[param2] = False
    
#     like=np.array(like)
#     like=like-np.amin(like)
    
#     return values1,values2,like

# def param_profile(m,param,min_val,max_val,size=50,iterate=5,subtract_min=True):
#     values=np.linspace(min_val,max_val,size)
#     like=[]
#     values_aux=[]
    
#     intial_values=np.array(m.values)
#     m.fixed[param] = True
    
#     for value in values:
# #         m.values=intial_values
#         m.values[param]=value
# #         m.scipy(method='SLSQP').migrad(iterate = iterate)
#         m.migrad(iterate = iterate)
#         if not m.valid:
#             print("Faield to converge in ",param," = ",value)
#         else:
#             like.append(m.fval)
#             values_aux.append(value)
#     m.fixed[param] = False
    
#     like=np.array(like)
#     if subtract_min==True:
#         like=like-np.amin(like)
    
#     return values_aux,like

# def param_profile_simplex(m,param,min_val,max_val,size=50):
#     values=np.linspace(min_val,max_val,size)
#     like=[]
    
#     intial_values=np.array(m.values)
#     m.fixed[param] = True
    
#     for value in values:
#         m.values=intial_values
#         m.values[param]=value
#         m.simplex.migrad()
#         if not m.valid:
#             print("Faield to converge in ",param," = ",value)
#         like.append(m.fval)
#     m.fixed[param] = False
    
#     like=np.array(like)
#     like=like-np.amin(like)
    
#     return values,like