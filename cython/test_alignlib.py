import alignlib_lite as alignlib

x = alignlib.py_makeAlignmentBlocks()

x.addDiagonal( 10, 100, 0 )

print x.getNumAligned(), x.getRowFrom(), x.getRowTo()


f = alignlib.py_AlignmentFormatBlat(x)
print str(f)
f.copy(x)
print str(f)


ff = alignlib.py_AlignmentFormatBlat("1\t100\t1001\t1101\t100,\t100,\t100,")
print str(ff)
ff.copy(x)
print str(x)
