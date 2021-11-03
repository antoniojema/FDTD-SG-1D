from GLOBAL import *

interp_old = eval(INTERP_OLD)
interp_new = eval(INTERP_NEW)

#####################
###   ADVANCE E   ###
#####################
def advanceE(E, H, CE, CH, CE_border, H_border, fine_range, n_E, sg, max_SG, from_E):
    E[sg][1:-1] = E[sg][1:-1] + CE[sg] * (H[sg][:-1] - H[sg][1:])
    n_E[sg] += 1

    if sg < max_SG:
        if (not BORDER_LTS) or (BORDER_LTS and n_E[sg] % BORDER_LTS_LOOP != 0):
            E[sg][ 0] = E[sg][ 0] + CE_border[sg] * (H_border[sg][ 0] -        H[sg][0])
            E[sg][-1] = E[sg][-1] + CE_border[sg] * (       H[sg][-1] - H_border[sg][1])

        if from_E:
            communicateToUpper(E, fine_range, sg)

    if sg > 0:
        if LTS:
            advanceE(E, H, CE, CH, CE_border, H_border, fine_range, n_E, sg-1, max_SG, from_E=True)
            advanceH(E, H, CE, CH, CE_border, H_border, fine_range, n_E, sg-1, max_SG)
        else:
            advanceE(E, H, CE, CH, CE_border, H_border, fine_range, n_E, sg-1, max_SG, from_E=True)


#####################
###   ADVANCE H   ###
######################
def advanceH(E, H, CE, CH, CE_border, H_border, fine_range, n_E, sg, max_SG):
    H[sg][:] = H[sg][:] + CH[sg] * (E[sg][:-1] - E[sg][1:])

    if sg > 0:
        if INTERPOLATING:
            communicateToLower(H, H_border, fine_range, sg, interpolate=True)
            pass
        else:
            communicateToLower(H, H_border, fine_range, sg, interpolate=False)
            pass

        if LTS:
            advanceE(E, H, CE, CH, CE_border, H_border, fine_range, n_E, sg-1, max_SG, from_E=False)
            advanceH(E, H, CE, CH, CE_border, H_border, fine_range, n_E, sg-1, max_SG)
        else:
            advanceH(E, H, CE, CH, CE_border, H_border, fine_range, n_E, sg-1, max_SG)

        if INTERPOLATING:
            communicateToLower(H, H_border, fine_range, sg, interpolate=False)
            pass
        else:
            #communicateToLower(H, H_border, fine_range, sg, interpolate=False)
            pass


###############
###   MUR   ###
###############
def advanceMur(E, E1, E2, CMur, max_SG):
    E[max_SG][ 0] = E1 + CMur * (E[max_SG][ 1] - E[max_SG][ 0])
    E[max_SG][-1] = E2 + CMur * (E[max_SG][-2] - E[max_SG][-1])

    E1 = E[max_SG][ 1]
    E2 = E[max_SG][-2]

    return E1, E2

################################
###   COMMUNICATE TO UPPER   ###
################################
def communicateToUpper(E, fine_range, sg):
    E[sg+1][fine_range[sg][0]] = E[sg][ 0]
    E[sg+1][fine_range[sg][1]] = E[sg][-1]

################################
###   COMMUNICATE TO LOWER   ###
################################
def communicateToLower(H, H_border, fine_range, sg, interpolate):
    if interpolate:
        H_border[sg-1][0] = interp_old * H_border[sg-1][0] + interp_new * H[sg][fine_range[sg-1][0]-1]
        H_border[sg-1][1] = interp_old * H_border[sg-1][1] + interp_new * H[sg][fine_range[sg-1][1]  ]
    else:
        H_border[sg-1][0] = H[sg][fine_range[sg-1][0]-1]
        H_border[sg-1][1] = H[sg][fine_range[sg-1][1]  ]
