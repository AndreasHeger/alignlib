import timeit
import alignlib

NUM_SAMPLES=1000
ALISIZE=2000


alignlib_vector = alignlib.makeAlignmentVector()
alignlib_vector.addDiagonal( 0, ALISIZE, 0)

python_vector = []
for x in xrange(ALISIZE): 
    python_vector.append(x) 

def pythonBuildVector():
    """build vector alignment in python."""
    vector = []
    for x in xrange(ALISIZE): 
        vector.append(x) 
     
def alignlibBuildVector():
    "Stupid test function"
    vector = alignlib.makeAlignmentVector()
    vector.addDiagonal( 0, ALISIZE, 0)

def pythonMapVector():
    "test speed of mapRowToCol"    
    for x in xrange(ALISIZE): 
        a = python_vector[x]
     
def alignlibMapVector():
    "test speed of mapRowToCol"
    for x in xrange(ALISIZE): 
        a = alignlib_vector.mapRowToCol(x)

def pythonMapVector():
    """build vector alignment in python."""
    for x in xrange(ALISIZE): 
        a = python_vector[x]
     
def alignlibCombineVector():
    "test combination of vectors"
    vector = alignlib.makeAlignmentVector()
    alignlib.combineAlignment( vector, alignlib_vector, alignlib_vector, alignlib.RR)

def pythonCombineVector():
    "test combination of vectors"
    x, y = 0, 0
    new = [0] * max(len(python_vector), len(python_vector))
    while x < len(python_vector) and y < len(python_vector):
        if x < y:
            x += 1
            continue
        elif y < x:
            y += 1
            continue
        mapped_x, mapped_y = python_vector[x], python_vector[y]
        new[mapped_x] = mapped_y
        x += 1
        y += 1

if __name__=='__main__':
    
    tests = ( ("pythonBuildVector", "alignlibBuildVector"), 
              ("pythonMapVector", "alignlibMapVector"),
              ("pythonCombineVector", "alignlibCombineVector" ) ) 
    
    for test1, test2 in tests:
        t1 = timeit.Timer("%s()" % test1, "from __main__ import %s" % test1)
        t2 = timeit.Timer("%s()" % test2, "from __main__ import %s" % test2)
        tt1 = t1.timeit( NUM_SAMPLES )
        tt2 = t2.timeit( NUM_SAMPLES )
        print "%s\t%s\t%5.2f\t%f\t%f" % (test2, test1, 100.0 * tt2 / tt1, tt2, tt1)
        
        
                 