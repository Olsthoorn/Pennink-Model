oc = [model, loading_package=False, budget_filerecord=None, concentration_filerecord=None, concentrationprintrecord=None, saverecord=None, printrecord=None, filename=None, pname=None, parent_file=None]


for p in oc:
    if '='  in p:
        p.split('=')
        print(p[0], '\t', p[1])