
#from GlobalsCorrector import *



def cprint12():
    global charly2, charly
    get_main_vars()
    print("cp",charly)
    print("cp",charly2)
    charly2="c"
    print("cp",charly,charly2)
    return None

def cprint3():
    global 3charly3
    get_main_vars()
    print("cp3",charly3)
    
    charly3=333
    print("cp3",charly3)
    return None



def get_main_vars():
    global charly2, charly,charly3
    from __main__ import charly2, charly,charly3
    print("inside get",charly,charly2)
    return None

print("GlobalsM loaded")