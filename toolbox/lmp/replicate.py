from data import data
from numpy import *
from itertools import combinations

# atom style must be full
# atoms have to be wrapped

replicate = [2,2,2]
inFile = 'cl.240'
outFile = 'rep.data'

d = data(inFile)

Lx_org = diff(d.headers["xlo xhi"])
Ly_org = diff(d.headers["ylo yhi"])
Lz_org = diff(d.headers["zlo zhi"])
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
    if irep == 0:   tmp[:,4] += (i*Lx)
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
  d.headers["atoms"] *= nrep
  d.headers["bonds"] *= nrep
  d.headers["angles"] *= nrep
  d.headers["dihedrals"] *= nrep
  d.headers["impropers"] *= nrep
  
  print d.headers["atoms"]

  if irep == 0:   d.headers["xlo xhi"] = (d.headers["xlo xhi"][0], float(d.headers["xlo xhi"][0]+(nrep*Lx)))
  elif irep == 1: d.headers["ylo yhi"] = (d.headers["ylo yhi"][0], float(d.headers["ylo yhi"][0]+(nrep*Ly)))
  else:           d.headers["zlo zhi"] = (d.headers["zlo zhi"][0], float(d.headers["zlo zhi"][0]+(nrep*Lz)))
    
  d.sections["Atoms"] = ['%d %d %d %f %f %f %f\n' % (a[0],a[1],a[2],a[3],a[4],a[5],a[6]) for a in atoms_new]
  d.sections["Bonds"] = ['%d %d %d %d\n' % (b[0],b[1],b[2],b[3]) for b in bonds_new]
  d.sections["Angles"] = ['%d %d %d %d %d\n' % (a[0],a[1],a[2],a[3],a[4]) for a in angle_new]
  d.sections["Dihedrals"] = ['%d %d %d %d %d %d\n' % (i[0],i[1],i[2],i[3],i[4],i[5]) for i in dihed_new]
  d.sections["Impropers"] = ['%d %d %d %d %d %d\n' % (i[0],i[1],i[2],i[3],i[4],i[5]) for i in impro_new]
  
  
# adjusting connectivity at periodic boundaries

atoms = array(d.get('Atoms'))
bonds = array(d.get('Bonds'))
angle = array(d.get('Angles'))
dihed = array(d.get('Dihedrals'))
impro = array(d.get('Impropers')) 

Natoms = d.headers["atoms"]
Nbonds = d.headers["bonds"]
Nangle = d.headers["angles"]
Ndihed = d.headers["dihedrals"]
Nimpro = d.headers["impropers"]

pbc = array([Lx,Ly,Lz])

for b in bonds:
  ib,t,i,j = map(int,b)
  if mod(ib,1000) == 0:print "bonds : ",ib,"/",Nbonds
  xi = atoms[atoms[:,0] == i,[4,5,6]]
  xj = atoms[atoms[:,0] == j,[4,5,6]]
  
  dij = linalg.norm(xj-xi)
  
  if any(dij > pbc/2.): 
    kmin = j - int(floor(j / float(Natoms_org))) * Natoms_org
    if kmin == 0: kmin = Natoms_org
    candidates = arange(kmin,Natoms+1,Natoms_org)   
    # finding closest atom to i from list of candidates
    dmin = linalg.norm(dij)
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
    
    
for a in angle:
  ia,t,a1,a2,a3 = map(int,a)
  if mod(ia,1000) == 0: print "angles: ",ia,"/",Nangle
  for i,j in combinations([a1,a2,a3],2):
    xi = atoms[atoms[:,0] == i,[4,5,6]]
    xj = atoms[atoms[:,0] == j,[4,5,6]]
    
    dij = linalg.norm(xj-xi)
    
    if any(dij > pbc/2.): 
      kmin = j - int(floor(j / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      # finding closest atom to i from list of candidates
      dmin = linalg.norm(dij)
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
      if j == a1: angle[angle[:,0] == ia,2] = kmin
      if j == a2: angle[angle[:,0] == ia,3] = kmin
      if j == a3: angle[angle[:,0] == ia,4] = kmin
     
      
for i in dihed:
  ii,t,i1,i2,i3,i4 = map(int,i)
  if mod(ii,1000) == 0: print "diheds: ",ii,"/",Ndihed
  
  if len(unique(i[2:])) < 4:
    print "WARNING: dihedral in a ring"
    
  
  for i,j in combinations([i1,i2,i3,i4],2):
    xi = atoms[atoms[:,0] == i,[4,5,6]]
    xj = atoms[atoms[:,0] == j,[4,5,6]]
    
    dij = linalg.norm(xj-xi)
    
    if any(dij > pbc/2.): 
      kmin = j - int(floor(j / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      # finding closest atom to i from list of candidates
      dmin = linalg.norm(dij)
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
      if j == i1: dihed[dihed[:,0] == ii,2] = kmin
      if j == i2: dihed[dihed[:,0] == ii,3] = kmin
      if j == i3: dihed[dihed[:,0] == ii,4] = kmin
      if j == i4: dihed[dihed[:,0] == ii,5] = kmin


for i in impro:
  ii,t,i1,i2,i3,i4 = map(int,i)
  if mod(ii,1000) == 0: print "impros: ",ii,"/",Nimpro
  
  if len(unique(i[2:])) < 4:
    print "WARNING: improper in a ring"
    
  for i,j in combinations([i1,i2,i3,i4],2):
    xi = atoms[atoms[:,0] == i,[4,5,6]]
    xj = atoms[atoms[:,0] == j,[4,5,6]]
    
    dij = linalg.norm(xj-xi)
    
    if any(dij > pbc/2.): 
      kmin = j - int(floor(j / float(Natoms_org))) * Natoms_org
      if kmin == 0: kmin = Natoms_org
      candidates = arange(kmin,Natoms+1,Natoms_org)   
      # finding closest atom to i from list of candidates
      dmin = linalg.norm(dij)
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
      if j == i1: impro[impro[:,0] == ii,2] = kmin
      if j == i2: impro[impro[:,0] == ii,3] = kmin
      if j == i3: impro[impro[:,0] == ii,4] = kmin
      if j == i4: impro[impro[:,0] == ii,5] = kmin

d.sections["Bonds"] = ['%d %d %d %d\n' % (b[0],b[1],b[2],b[3]) for b in bonds]   
d.sections["Angles"] = ['%d %d %d %d %d\n' % (a[0],a[1],a[2],a[3],a[4]) for a in angle]  
d.sections["Dihedrals"] = ['%d %d %d %d %d %d\n' % (a[0],a[1],a[2],a[3],a[4],a[5]) for a in dihed]  
d.sections["Impropers"] = ['%d %d %d %d %d %d\n' % (a[0],a[1],a[2],a[3],a[4],a[5]) for a in impro]   

# remove velocities
d.sections["Velocities"] = []

d.write(outFile)
