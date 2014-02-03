import alignlib_lite as alignlib

x = alignlib.py_makeAlignmentBlocks()

x.addDiagonal( 10, 100, 0 )

print x.getNumAligned(), x.getRowFrom(), x.getRowTo()


f = alignlib.py_AlignmentFormatBlat( x )

print str(f)

f.copy( x )

print str(f)
