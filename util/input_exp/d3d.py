def shotinfo(shot):
    # LTO 1: after 210 counter, 2: after 150 off-axis 3: after 210 off-axis
    if shot > 180000: LTO = 3 # <== need to update
    elif shot > 143703: LTO = 2 
    elif shot > 125042: LTO = 1
    else: LTO = 0
    return LTO

