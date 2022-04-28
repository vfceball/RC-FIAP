import openseespy.opensees as op
def Exci_pattern(dt, outFile, sf):
    # Uniform EXCITATION: acceleration input
    op.timeSeries('Path', 2, '-dt', dt, '-filePath', outFile, '-factor', sf)
    op.pattern('UniformExcitation', 400, 1, '-accel', 2)
