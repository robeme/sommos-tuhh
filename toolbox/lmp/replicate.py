from data import data
from numpy import *

# atom style must be full
# atoms have to be wrapped

replicate = [2,2,2]
inFile = 'cl.240'
outFile = 'rep.data'

d = data(inFile)

Lx_org = diff(d.headers["xlo xhi"])
Ly_org = diff(d.headers["ylo yhi"])
Lz_org = diff(d.headers["zlo zhi"])
pbc_org = array([Lx_org,Ly_org,Lz_org])
Natoms_org = d.headers["atoms"]

for irep,nrep in enumerate(replicate):

  Natoms = d.headers["atoms"]
  Nbonds = d.headers["bonds"]
  Nangle = d.headers["angles"]
  Ndihed = d.headers["dihedrals"]
  Nimpro = d.headers["impropers"]
  
  Nmols = max(d.get('Atoms',2))

  Lx = diff(d.headers["xlo xhi"])
  Ly = diff(d.headers["ylo yhi"])
  Lz = diff(d.headers["zlo zhi"])

  atoms = array(d.get('Atoms'))
  bonds = array(d.get('Bonds'))
  angle = array(d.get('Angles'))
  dihed = array(d.get('Dihedrals'))
  impro = array(d.get('Impropers')) 

  atoms_new = array([]).reshape(-1,atoms.shape[1])
  bonds_new = array([]).reshape(-1,bonds.shape[1])
  angle_new = array([]).reshape(-1,angle.shape[1])
  dihed_new = array([]).reshape(-1,dihed.shape[1])
  impro_new = array([]).reshape(-1,impro.shape[1])

  for i in range(nrep):
    tmp = copy(atoms)
    # change atomic ids
    tmp[:,0] += (i*Natoms)
    # and keep track of the molecule ids
    # I do not take care of broken bonds maybe I should
    # recreate the molecules afterwards again from bond
    # topology
    tmp[:,1] += (i*Nmols)
    # add offset in x
    if   irep == 0: tmp[:,4] += (i*Lx)
    elif irep == 1: tmp[:,5] += (i*Ly)
    else:           tmp[:,6] += (i*Lz)
    atoms_new = vstack([atoms_new,tmp])
    
    tmp = copy(bonds)
    # increase numbering of bonds
    tmp[:,0] += (i*Nbonds)
    # change atomic ids
    tmp[:,2] += (i*Natoms)
    tmp[:,3] += (i*Natoms)
    bonds_new = vstack([bonds_new,tmp])
    
    tmp = copy(angle)
    # increase numbering of angles
    tmp[:,0] += (i*Nangle)
    # change atomic ids
    tmp[:,2] += (i*Natoms)
    tmp[:,3] += (i*Natoms)
    tmp[:,4] += (i*Natoms)
    angle_new = vstack([angle_new,tmp])
    
    tmp = copy(dihed)
    # increase numbering of dihedrals
    tmp[:,0] += (i*Ndihed)
    # change atomic ids
    tmp[:,2] += (i*Natoms)
    tmp[:,3] += (i*Natoms)
    tmp[:,4] += (i*Natoms)
    tmp[:,5] += (i*Natoms)
    dihed_new = vstack([dihed_new,tmp])
    
    tmp = copy(impro)
    # increase numbering of impropers
    tmp[:,0] += (i*Nimpro)
    # change atomic ids
    tmp[:,2] += (i*Natoms)
    tmp[:,3] += (i*Natoms)
    tmp[:,4] += (i*Natoms)
    tmp[:,5] += (i*Natoms)
    impro_new = vstack([impro_new,tmp])
    
  # reset headers and counters
  d.headers["atoms"]     *= nrep
  d.headers["bonds"]     *= nrep
  d.headers["angles"]    *= nrep
  d.headers["dihedrals"] *= nrep
  d.headers["impropers"] *= nrep

  if   irep == 0: d.headers["xlo xhi"] = (d.headers["xlo xhi"][0], float(d.headers["xlo xhi"][0]+(nrep*Lx)))
  elif irep == 1: d.headers["ylo yhi"] = (d.headers["ylo yhi"][0], float(d.headers["ylo yhi"][0]+(nrep*Ly)))
  else:           d.headers["zlo zhi"] = (d.headers["zlo zhi"][0], float(d.headers["zlo zhi"][0]+(nrep*Lz)))
    
  d.sections["Atoms"]     = ['%d %d %d %.12f %f %f %f\n' % (a[0],a[1],a[2],a[3],a[4],a[5],a[6]) for a in atoms_new]
  d.sections["Bonds"]     = ['%d %d %d %d\n'             % (b[0],b[1],b[2],b[3])                for b in bonds_new]
  d.sections["Angles"]    = ['%d %d %d %d %d\n'          % (a[0],a[1],a[2],a[3],a[4])           for a in angle_new]
  d.sections["Dihedrals"] = ['%d %d %d %d %d %d\n'       % (i[0],i[1],i[2],i[3],i[4],i[5])      for i in dihed_new]
  d.sections["Impropers"] = ['%d %d %d %d %d %d\n'       % (i[0],i[1],i[2],i[3],i[4],i[5])      for i in impro_new]
  
  
# adjusting connectivity at periodic boundaries

atoms = array(d.get('Atoms'))
bonds = array(d.get('Bonds'))
angle = array(d.get('Angles'))
dihed = array(d.get('Dihedrals'))
impro = array(d.get('Impropers')) 

Lx = diff(d.headers["xlo xhi"])
Ly = diff(d.headers["ylo yhi"])
Lz = diff(d.headers["zlo zhi"])

Natoms = d.headers["atoms"]
Nbonds = d.headers["bonds"]
Nangle = d.headers["angles"]
Ndihed = d.headers["dihedrals"]
Nimpro = d.headers["impropers"]

pbc = array([Lx,Ly,Lz])

overstretched = []
print "correcting bonds ..."
for b in bonds:
  ib,t,i,j = map(int,b)
  if mod(ib,1000) == 0:print "  ",ib,"/",Nbonds
  xi = atoms[atoms[:,0] == i,[4,5,6]]
  xj = atoms[atoms[:,0] == j,[4,5,6]]
  
  xij = xj-xi
  
  if any(abs(xij) > pbc_org/2.):
    # found overstretched bond and create list of possible other neighbors
    kmin = j - int(floor(j / float(Natoms_org))) * Natoms_org
    if kmin == 0: kmin = Natoms_org
    candidates = arange(kmin,Natoms+1,Natoms_org)   
    #finding closest atom to i from list of candidates
    dmin = linalg.norm(xij)
    kmin = j
    for k in candidates:
      xk = atoms[atoms[:,0] == k,[4,5,6]]
      xk = xk - floor(xk / pbc) * pbc
      dx = xk - xi
      dx = linalg.norm(dx - around(dx / pbc) * pbc)
      if dx < dmin: 
        dmin = dx
        kmin = k
    # reset bond connectivity to closest replicated atom
    bonds[bonds[:,0] == ib,3] = kmin
    # store atoms in each overstretched bond for later 
    # use in correcting angles, dihedrals and impropers
    overstretched.append([i,j])

# find overstretched bonds in angles, diheds and impros
print "correcting angles, diheds, impros ..."
for ib,ob in enumerate(overstretched):
  if mod(ib,100) == 0:print "  ",ib,"/",len(overstretched)
  # loop over angles ids which have an overstretched bond
  for idx in angle[(any(angle[:,[2,3,4]] == ob[0],axis=1) & any(angle[:,[2,3,4]] == ob[1],axis=1)),0]:
    a1,a2,a3 = angle[angle[:,0]==idx,[2,3,4]]
    
    xi = atoms[atoms[:,0] == a1,[4,5,6]]
    xj = atoms[atoms[:,0] == a2,[4,5,6]]
  
    xij = xj-xi
    
    if any(abs(xij) > pbc_org/2.):
      kmin = a2 - int(floor(a2 / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      dmin = linalg.norm(xij)
      kmin = a2
      for k in candidates:
        xk = atoms[atoms[:,0] == k,[4,5,6]]
        xk = xk - floor(xk / pbc) * pbc
        dx = xk - xi
        dx = linalg.norm(dx - around(dx / pbc) * pbc)
        if dx < dmin: 
          dmin = dx
          kmin = k
      a2 = kmin
    
    xi = atoms[atoms[:,0] == a2,[4,5,6]]
    xj = atoms[atoms[:,0] == a3,[4,5,6]]
  
    xij = xj-xi
    
    if any(abs(xij) > pbc_org/2.):
      kmin = a3 - int(floor(a3 / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      dmin = linalg.norm(xij)
      kmin = a3
      for k in candidates:
        xk = atoms[atoms[:,0] == k,[4,5,6]]
        xk = xk - floor(xk / pbc) * pbc
        dx = xk - xi
        dx = linalg.norm(dx - around(dx / pbc) * pbc)
        if dx < dmin: 
          dmin = dx
          kmin = k
      a3 = kmin
      
    # reset overstretched angle
    angle[angle[:,0]==idx,2] = a1
    angle[angle[:,0]==idx,3] = a2 
    angle[angle[:,0]==idx,4] = a3

  # loop over diheds ids which have an overstretched bond
  for idx in dihed[(any(dihed[:,[2,3,4,5]] == ob[0],axis=1) & any(dihed[:,[2,3,4,5]] == ob[1],axis=1)),0]:
    a1,a2,a3,a4 = dihed[dihed[:,0]==idx,[2,3,4,5]]
    
    xi = atoms[atoms[:,0] == a1,[4,5,6]]
    xj = atoms[atoms[:,0] == a2,[4,5,6]]
  
    xij = xj-xi
    
    if any(abs(xij) > pbc_org/2.):
      kmin = a2 - int(floor(a2 / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      dmin = linalg.norm(xij)
      kmin = a2
      for k in candidates:
        xk = atoms[atoms[:,0] == k,[4,5,6]]
        xk = xk - floor(xk / pbc) * pbc
        dx = xk - xi
        dx = linalg.norm(dx - around(dx / pbc) * pbc)
        if dx < dmin: 
          dmin = dx
          kmin = k
      a2 = kmin
    
    xi = atoms[atoms[:,0] == a2,[4,5,6]]
    xj = atoms[atoms[:,0] == a3,[4,5,6]]
  
    xij = xj-xi
    
    if any(abs(xij) > pbc_org/2.):
      kmin = a3 - int(floor(a3 / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      dmin = linalg.norm(xij)
      kmin = a3
      for k in candidates:
        xk = atoms[atoms[:,0] == k,[4,5,6]]
        xk = xk - floor(xk / pbc) * pbc
        dx = xk - xi
        dx = linalg.norm(dx - around(dx / pbc) * pbc)
        if dx < dmin: 
          dmin = dx
          kmin = k
      a3 = kmin
      
    xi = atoms[atoms[:,0] == a3,[4,5,6]]
    xj = atoms[atoms[:,0] == a4,[4,5,6]]
  
    xij = xj-xi
    
    if any(abs(xij) > pbc_org/2.):
      kmin = a4 - int(floor(a4 / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      dmin = linalg.norm(xij)
      kmin = a4
      for k in candidates:
        xk = atoms[atoms[:,0] == k,[4,5,6]]
        xk = xk - floor(xk / pbc) * pbc
        dx = xk - xi
        dx = linalg.norm(dx - around(dx / pbc) * pbc)
        if dx < dmin: 
          dmin = dx
          kmin = k
      a4 = kmin
      
    # reset overstretched dihedral
    dihed[dihed[:,0]==idx,2] = a1
    dihed[dihed[:,0]==idx,3] = a2 
    dihed[dihed[:,0]==idx,4] = a3
    dihed[dihed[:,0]==idx,5] = a4     
    
  # loop over impro ids which have an overstretched bond
  # the order of the atoms depends in lammps on the improper
  # style used and actually for an improper there must not 
  # necessarily a bond between i-j,j-k and k-l defined. However,
  # typically a bond exist between i-j,i-k and i-l is defined
  # in an improper.
  for idx in impro[(any(impro[:,[2,3,4,5]] == ob[0],axis=1) & any(impro[:,[2,3,4,5]] == ob[1],axis=1)),0]:
    a1,a2,a3,a4 = impro[impro[:,0]==idx,[2,3,4,5]]
    
    xi = atoms[atoms[:,0] == a1,[4,5,6]]
    xj = atoms[atoms[:,0] == a2,[4,5,6]]
  
    xij = xj-xi
    
    if any(abs(xij) > pbc_org/2.):
      kmin = a2 - int(floor(a2 / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      dmin = linalg.norm(xij)
      kmin = a2
      for k in candidates:
        xk = atoms[atoms[:,0] == k,[4,5,6]]
        xk = xk - floor(xk / pbc) * pbc
        dx = xk - xi
        dx = linalg.norm(dx - around(dx / pbc) * pbc)
        if dx < dmin: 
          dmin = dx
          kmin = k
      a2 = kmin
    
    xi = atoms[atoms[:,0] == a1,[4,5,6]]
    xj = atoms[atoms[:,0] == a3,[4,5,6]]
  
    xij = xj-xi
    
    if any(abs(xij) > pbc_org/2.):
      kmin = a3 - int(floor(a3 / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      dmin = linalg.norm(xij)
      kmin = a3
      for k in candidates:
        xk = atoms[atoms[:,0] == k,[4,5,6]]
        xk = xk - floor(xk / pbc) * pbc
        dx = xk - xi
        dx = linalg.norm(dx - around(dx / pbc) * pbc)
        if dx < dmin: 
          dmin = dx
          kmin = k
      a3 = kmin
      
    xi = atoms[atoms[:,0] == a1,[4,5,6]]
    xj = atoms[atoms[:,0] == a4,[4,5,6]]
  
    xij = xj-xi
    
    if any(abs(xij) > pbc_org/2.):
      kmin = a4 - int(floor(a4 / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      dmin = linalg.norm(xij)
      kmin = a4
      for k in candidates:
        xk = atoms[atoms[:,0] == k,[4,5,6]]
        xk = xk - floor(xk / pbc) * pbc
        dx = xk - xi
        dx = linalg.norm(dx - around(dx / pbc) * pbc)
        if dx < dmin: 
          dmin = dx
          kmin = k
      a4 = kmin
      
    # reset overstretched dihedral
    impro[impro[:,0]==idx,2] = a1
    impro[impro[:,0]==idx,3] = a2 
    impro[impro[:,0]==idx,4] = a3
    impro[impro[:,0]==idx,5] = a4 

d.sections["Bonds"] = ['%d %d %d %d\n' % (b[0],b[1],b[2],b[3]) for b in bonds]   
d.sections["Angles"] = ['%d %d %d %d %d\n' % (a[0],a[1],a[2],a[3],a[4]) for a in angle]  
d.sections["Dihedrals"] = ['%d %d %d %d %d %d\n' % (a[0],a[1],a[2],a[3],a[4],a[5]) for a in dihed]  
d.sections["Impropers"] = ['%d %d %d %d %d %d\n' % (a[0],a[1],a[2],a[3],a[4],a[5]) for a in impro]   

# remove velocities
d.delete("Velocities")

d.write(outFile)
