from Eout import calculate_Eout
from kinematics import kinematics

def e_after_JENSA(reaction, theta, physicalTgt, stripfoil, dedx_lib):
    ion = reaction.beam
    e_ion_after_physicalTgt_list = []
    e_ion_after_stripfoil_list = []
    if physicalTgt.status:
        if reaction.status:
            rxn_distance = physicalTgt.thickness / 2
            e_ion_at_rxn = calculate_Eout(ion, reaction.eb * 1000, physicalTgt.material, rxn_distance, dedx_lib)
            ion = reaction.recoil
            reacKin = kinematics(reaction.m1, reaction.m2, reaction.m3,
                                         reaction.m4, reaction.exc3, e_ion_at_rxn / 1000, theta)  # energy in GeV
            reacKin.pop()
            for e_ion_after_rxn in reacKin:
                e_ion_after_physicalTgt_list.append(calculate_Eout(ion, e_ion_after_rxn * 1000, physicalTgt.material, physicalTgt.thickness - rxn_distance, dedx_lib))
        else:
            e_ion_after_physicalTgt_list.append(calculate_Eout(ion, reaction.eb * 1000, physicalTgt.material, physicalTgt.thickness, dedx_lib))
    else:
        e_ion_after_physicalTgt_list.append(reaction.eb * 1000)

    for e_ion_after_physicalTgt in e_ion_after_physicalTgt_list:
        if stripfoil.status:
            e_ion_after_stripfoil_list.append(calculate_Eout(ion, e_ion_after_physicalTgt, stripfoil.material, stripfoil.thickness, dedx_lib))
        else:
            e_ion_after_stripfoil_list.append(e_ion_after_physicalTgt)

    return e_ion_after_stripfoil_list


def get_pid(ion, e_ion_after_stripfoil, mcp, ic, dssd, dedx_lib):

    if mcp.status:
        e_ion_after_mcp = calculate_Eout(ion, e_ion_after_stripfoil, mcp.material, mcp.thickness, dedx_lib)
    else:
        e_ion_after_mcp = e_ion_after_stripfoil

    e_ion_after_ic_window = calculate_Eout(ion, e_ion_after_mcp, ic.window_material, ic.window_thickness, dedx_lib)
    if ic.gas_status:
        e_ion_after_ic_dl1 = calculate_Eout(ion, e_ion_after_ic_window, ic.gas_material, ic.dl1_thickness, dedx_lib)
        e_ion_after_ic_de = calculate_Eout(ion, e_ion_after_ic_dl1, ic.gas_material, ic.de_thickness, dedx_lib)
        e_ion_after_ic_dl2 = calculate_Eout(ion, e_ion_after_ic_de, ic.gas_material, ic.dl2_thickness, dedx_lib)
        delE = e_ion_after_ic_dl1 - e_ion_after_ic_de
    else:
        e_ion_after_ic_dl2 = e_ion_after_ic_window
        delE = 0

    e_ion_after_dssd_dl = calculate_Eout(ion, e_ion_after_ic_dl2, dssd.dl_material, dssd.dl_thickness, dedx_lib)
    e_ion_after_dssd = calculate_Eout(ion, e_ion_after_dssd_dl, dssd.det_material, dssd.det_thickness, dedx_lib)
    E = e_ion_after_dssd_dl - e_ion_after_dssd

    pid = [delE, E]

    # print(e_ion_after_physicalTgt, e_ion_after_mcp, e_ion_after_ic_window, e_ion_after_ic_dl1, e_ion_after_ic_de, e_ion_after_ic_dl2, e_ion_after_dssd_dl, e_ion_after_dssd)

    return pid