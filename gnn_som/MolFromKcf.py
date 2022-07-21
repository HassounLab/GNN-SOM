from rdkit import Chem
from rdkit.Chem.rdmolops import FastFindRings
from rdkit.Geometry.rdGeometry import Point3D

def MolFromKcfContents(contents):
    def parseChargeOrRadical(field):
        assert field[0] == '#'
        count = field[1:-1]
        if len(count) == 0:
            count = 1 # implicit 1
        else:
            count = int(count)
        if field[-1] == '-':
            return ('charge', -count)
        elif field[-1] == '+':
            return ('charge', count)
        elif field[-1] == '^':
            return ('radical', count)
        else:
            assert False
        return charge

    '''
    assert parseChargeOrRadical('#+') == ('charge', 1)
    assert parseChargeOrRadical('#-') == ('charge', -1)
    assert parseChargeOrRadical('#3+') == ('charge', 3)
    assert parseChargeOrRadical('#^') == ('radical', 1)
    '''
    
    bondTypes = [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]
    mol = Chem.RWMol()
    conformer = None
    atomIndices = {} # KCF atom index: RdKit index
    
    atoms = []
    bonds = []
    readingAtoms = False
    readingBonds = False
    for line in contents.splitlines():
        label = line[0:12].strip()
        content = line[12:].strip()
        
        if label == 'COMPOUND':
            mol.SetProp('_Name', content)
        
        if label == 'ATOM':
            atomCount = int(content)
            conformer = Chem.Conformer(atomCount)
            readingAtoms = True
            continue
        elif label != '':
            readingAtoms = False
        
        if label == 'BOND':
            readingBonds = True
            continue # skip bond count
        elif label != '':
            readingBonds = False
        
        if readingAtoms:
            fields = content.split()
            index, type, element, x, y = fields[0:5]
            index = int(index)
            x, y = float(x), float(y)
            tag = parseChargeOrRadical(fields[5]) if len(fields) > 5 else None
            
            if element == 'H+': # C00080
                element = 'H'
                assert tag is None # make sure we don't overwrite it
                tag = ('charge', 1)
            elif element == 'OH':
                element = 'O'
                assert tag is None
                tag = ('charge', -1)
            elif element in ['R', 'R#', 'X']:
                element = '*'
            
            atomIndices[index] = mol.AddAtom(Chem.Atom(element))
            if tag is not None:
                if tag[0] == 'charge':
                    mol.GetAtomWithIdx(atomIndices[index]).SetFormalCharge(tag[1])
                elif tag[0] == 'radical':
                    mol.GetAtomWithIdx(atomIndices[index]).SetNumRadicalElectrons(tag[1])
                else:
                    assert False
            mol.GetAtomWithIdx(atomIndices[index]).SetProp('kcfIndex', str(index))
            mol.GetAtomWithIdx(atomIndices[index]).SetProp('kcfType', type)
            conformer.SetAtomPosition(atomIndices[index], Point3D(x, y, 0))
        
        elif readingBonds:
            fields = content.split()
            index, atom1, atom2, order = fields[0:4]
            atom1, atom2, order = int(atom1), int(atom2), int(order)

            bondIndex = mol.AddBond(atomIndices[atom1], atomIndices[atom2], bondTypes[order - 1])
            bond = mol.GetBondWithIdx(bondIndex - 1) # are indices returned by AddBond off by one?
            if len(fields) == 4:
                bond.SetBondDir(Chem.BondDir.NONE)
            elif fields[4] == '#Up':
                bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif fields[4] == '#Down':
                bond.SetBondDir(Chem.BondDir.BEGINDASH)
            elif fields[4] == '#Either' and order == 1:
                bond.SetBondDir(Chem.BondDir.UNKNOWN)
            elif fields[4] == '#Either' and order == 2:
                # special case for "either" double bonds, see
                # https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/FileParsers/MolFileParser.cpp#L1483
                bond.SetBondDir(Chem.BondDir.EITHERDOUBLE)
                bond.SetStereo(Chem.BondStereo.STEREOANY)
            else:
                assert False
    
    # Setting atom position using a conformer:
    # https://sourceforge.net/p/rdkit/mailman/message/36474923/
    mol.AddConformer(conformer, assignId=True)
    
    # Needed for ExactMolWt(), CalcMolFormula(), etc.
    # https://loudspeaker.sakura.ne.jp/devblog/2019/01/22/way-rescue-compounds-not-loaded-smiles/
    mol.UpdatePropertyCache(strict=False)
    FastFindRings(mol)
    
    # Needed to correctly orient wedge/dash bonds
    Chem.AssignChiralTypesFromBondDirs(mol)
    
    # Resolve ambigous KEGG atom types for oxygens
    for atom in mol.GetAtoms():
        if atom.GetProp('kcfType') not in ['O6a', 'O7a', 'O7x', 'O1c']:
            continue
        
        bondCounts = {t: 0 for t in bondTypes}
        for bond in atom.GetBonds():
            bondCounts[bond.GetBondType()] += 1
        
        assert bondCounts[Chem.BondType.TRIPLE] == 0 # triple bonds are not possible for oxygen
        
        if bondCounts[Chem.BondType.DOUBLE] == 1: # attached via a double bond
            assert bondCounts[Chem.BondType.SINGLE] == 0 # single bonds not possible
            if atom.GetProp('kcfType') == 'O1c': # labeled P-OH but actually has a double bond (see C00013)
                atom.SetProp('kcfType', 'O3b') # P=O
            else: # O6a, O7a, O7x
                atom.SetProp('kcfType', atom.GetProp('kcfType') + '2')
        elif bondCounts[Chem.BondType.SINGLE] >= 1: # attached via a single bond
            atom.SetProp('kcfType', atom.GetProp('kcfType') + '1')
    return mol

def MolFromKcfFile(filename):
    with open(filename, 'r') as f:
        return MolFromKcfContents(f.read())
