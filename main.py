import matplotlib.pyplot as plt
from get_pid import *
from kinematics import *
import os
from get_mass import get_m
from setup_classes import *
import numpy as np


if __name__ == '__main__':

#======== Input = Instrument setup ========#
    physicalTgt = PhysicalTarget(1, 'gas', '240psi', 2500, 'helium')
    # status, state (gas/solid), pressure (0 for solid target, psi for physicalTgt, torr for gas cell), \
    # thickness (um), material
    # 2.5mm 150 = (8.48e-5g/cm3 3.19at/cm2) 160 = (9.1e-5g/cm3 3.42at/cm2)

    mcp = MCPs(0, 0.5, 'mylar')
    # status, thickness (um), material

    stripfoil = StripFoil(1, 0.2, 'carbon')
    # status, thickness (um), material
    # 0.2 um = 45 ug/cm2

    ic = IC(3, 'mylar', 1, '50T', 12700, 50800, 3200, 'isobutane')
    # window thickness (um), window material, gas status, gas pressure,
    # dl1 thickness (um), de thickness, dl2 thickness, gas material

    dssd = DSSD(0.5, 'aluminium', 300, 'silicon')
    # dl thickness, dl material (um), det thickness, det material

#======== Input = Reaction setup ========#
    target = 'He4' # particle1
    beam = 'Kr86' # particle2
    recoil = 'Sr88' # particle3. Traditional kinematics would call it ejectile. But in SECAR this is recoil
    ejectile = 'n2' # particle4. n0 for gamma.
    beam_q = 27
    recoil_q = 27
    beam_energy = 2.5*86 # MeV
    dedx_lib = 'SRIM' # this option requires user to calculate stopping power tables for ion in each material
    # dedx_lib = 'VICAR'

#======== Input = Plotting Parameters ========#
    n_rec = 100 # number of recoil particles generated in Monte-Carlo (theta_com and theta_lab)
    leaky_b_range = [0.8, 1.0] # Calculate leaky beam PID for beam energy fraction in this range.
    leaky_b_step = 0.01
    E_range = [70, 130] # DSSD energy axis range in MeV
    delE_range = [20, 100] # IC_dE energy axis range in MeV
    brho_accept = 0.0154 # +/- fraction of the recoil Brho value
    vel_accept = 0.0154 # +/- fraction of the recoil velocity value

# ======== End of user input =================#

# ======== Calculations ======================#
    # calculate masses of ions in reaction in GeV/c2
    m_list = [get_m(target), get_m(beam), get_m(recoil), get_m(ejectile)] # micro-amu
    for i in range(0, 4): m_list[i] = m_list[i]*931.494*1e-9 # GeV/c2

    rxn = Reaction(1, beam, recoil, m_list[0], m_list[1], m_list[2], m_list[3], 0.0, beam_energy*1e-3)
    # status, beam, recoil, m1 (GeV/c2), m2, m3, m4, exc3, eb (GeV)

    # compile the VICAR stopping power c++ code.
    if dedx_lib == 'VICAR': os.system('cd ./VICAR_dedx/ && g++ stp.cxx && cd ../')

    # Calculate and print theta_max in lab frame
    tmax = kinematics(rxn.m1, rxn.m2, rxn.m3, rxn.m4, rxn.exc3, rxn.eb, 0.000)[2]
    print ('theta_max = {:.2f} mrad'.format(tmax*1000))

    # Arrays for reaction kinematics plots
    theta_com = []
    theta_lab = []
    theta_lab2 = []
    energy_lab = []
    for i in range(0, n_rec, 1):
        # Generate random COM theta for an isotropic distribution in COM
        theta_com.append(generate_random_theta_com())
        # COM to lab conversion
        theta_lab.append(cm2lab_theta(theta_com[i], m_list[0], m_list[1], m_list[2], m_list[3], 0.0, beam_energy*1e-3))
        # Calculate energy solutions for recoil with theta = theta_lab
        e_kin = kinematics(rxn.m1, rxn.m2, rxn.m3, rxn.m4, rxn.exc3, rxn.eb, theta_lab[i])
        for j in range(0,2,1):
            theta_lab2.append(theta_lab[i])
            energy_lab.append(e_kin[j])

    # Recoil arrays
    delE_recoil = []
    E_recoil = []
    brho_recoil = []
    vel_recoil = []
    delE_recoil_brho = []
    E_recoil_brho = []
    delE_recoil_vel = []
    E_recoil_vel = []
    delE_recoil_brho_vel = []
    E_recoil_brho_vel = []

    # Leaky beam arrays
    delE_beam = []
    E_beam = []
    delE_beam_vel = []
    E_beam_vel = []
    delE_beam_delQminus2 = []
    E_beam_delQminus2 = []
    delE_beam_delQminus1 = []
    E_beam_delQminus1 = []
    delE_beam_delQ0 = []
    E_beam_delQ0 = []
    delE_beam_delQplus1 = []
    E_beam_delQplus1 = []
    delE_beam_delQplus2 = []
    E_beam_delQplus2 = []
    delE_beam_delQminus2_vel = []
    E_beam_delQminus2_vel = []
    delE_beam_delQminus1_vel = []
    E_beam_delQminus1_vel = []
    delE_beam_delQ0_vel = []
    E_beam_delQ0_vel = []
    delE_beam_delQplus1_vel = []
    E_beam_delQplus1_vel = []
    delE_beam_delQplus2_vel = []
    E_beam_delQplus2_vel = []

    # Set reaction and physicalTgt status ON.
    rxn.status = 1
    physicalTgt.status = 1

    # Filling recoil particle ID arrays
    for theta in theta_lab:
        # Calculate recoil's energy after reaction and (physicalTgt + stripfoil) thicknesses
        e_into_SECAR = e_after_JENSA(rxn, theta, physicalTgt, stripfoil, dedx_lib)
        for energy in e_into_SECAR:
            # Calculate brho of each recoil entering SECAR
            brho_recoil.append(get_brho(energy, rxn.m3, recoil_q))
            # Calculate velocity of each recoil entering SECAR
            vel_recoil.append(get_vel(energy, rxn.m3))
            # Calculate energy deposited by recoil in IC_dE and DSSD
            pid = get_pid(rxn.recoil, energy, mcp, ic, dssd, dedx_lib)
            delE_recoil.append(pid[0])
            E_recoil.append(pid[1])

    # Calculate mean brho and mean velocity
    brho_recoil_avg = sum(brho_recoil)/len(brho_recoil)
    vel_recoil_avg = sum(vel_recoil)/len(vel_recoil)

    for i in range(0,len(brho_recoil),1):
        if brho_recoil[i] < 1.015*brho_recoil_avg and brho_recoil[i] > 0.985*brho_recoil_avg:
            delE_recoil_brho.append(delE_recoil[i])
            E_recoil_brho.append(E_recoil[i])
            if vel_recoil[i] < 1.015 * vel_recoil_avg and vel_recoil[i] > 0.985 * vel_recoil_avg:
                delE_recoil_brho_vel.append(delE_recoil[i])
                E_recoil_brho_vel.append(E_recoil[i])
        if vel_recoil[i] < 1.015*vel_recoil_avg and vel_recoil[i] > 0.985*vel_recoil_avg:
            delE_recoil_vel.append(delE_recoil[i])
            E_recoil_vel.append(E_recoil[i])

    print('Average recoil brho and vel = {} and {}'.format(brho_recoil_avg, vel_recoil_avg))

    print('For an isotropic distribution of recoils in center of mass:')
    brho_acceptance = len(E_recoil_brho)/len(E_recoil)
    err_brho_acceptance = brho_acceptance*((1/len(E_recoil_brho))+(1/len(E_recoil)))**0.5
    print('\t Brho acceptance (+/- {:.1f}%) = {:.2f} +/- {:.2f} %'.format(100*brho_accept, 100*brho_acceptance, 100*err_brho_acceptance))

    vel_acceptance = len(E_recoil_vel) / len(E_recoil)
    err_vel_acceptance = vel_acceptance * ((1 / len(E_recoil_vel)) + (1 / len(E_recoil))) ** 0.5
    print('\t Velocity acceptance (+/- {:.1f}%) = {:.2f} +/- {:.2f} %'.format(100*vel_accept, 100 * vel_acceptance, 100 * err_vel_acceptance))

    brho_vel_acceptance = len(E_recoil_brho_vel) / len(E_recoil)
    err_brho_vel_acceptance = brho_vel_acceptance * ((1 / len(E_recoil_brho_vel)) + (1 / len(E_recoil))) ** 0.5
    print('\t Brho + velocity acceptance (+/- {:.1f}%) = {:.2f} +/- {:.2f} %'.format(100 * brho_accept, 100 * brho_vel_acceptance,
                                                                          100 * err_brho_vel_acceptance))

    # Turn OFF the reaction
    rxn.status = 0
    # Get leaky beam PID in energy range set by user in step sizes set by user
    eb_orig = rxn.eb
    for i in np.arange(leaky_b_range[0], leaky_b_range[1], leaky_b_step):
        rxn.eb = eb_orig * i
        # Calculate energy loss in physicalTgt and stripfoil.
        e_into_SECAR = e_after_JENSA(rxn, 0, physicalTgt, stripfoil, dedx_lib)[0]
        # Calculate energies deposited in IC_dE and DSSD
        pid = get_pid(rxn.beam, e_into_SECAR, mcp, ic, dssd, dedx_lib)
        delE_beam.append(pid[0])
        E_beam.append(pid[1])
        # calculate particle's velocity
        vel_beam = get_vel(e_into_SECAR, rxn.m2)
        # Check if the velocity matches the acceptance criteria set by user
        if vel_beam<(1+vel_accept)*vel_recoil_avg and vel_beam>(1-vel_accept)*vel_recoil_avg:
            delE_beam_vel.append(pid[0])
            E_beam_vel.append(pid[1])
        # Check brho of different charge states and compare with acceptance criteria
        for i in range(-2,3,1):
            brho_beam = get_brho(e_into_SECAR, rxn.m2, beam_q+i)
            if brho_beam < (1+brho_accept)*brho_recoil_avg and brho_beam > (1-brho_accept)*brho_recoil_avg:
                if i == -2:
                    delE_beam_delQminus2.append(pid[0])
                    E_beam_delQminus2.append(pid[1])
                elif i == -1:
                    delE_beam_delQminus1.append(pid[0])
                    E_beam_delQminus1.append(pid[1])
                elif i == 0:
                    delE_beam_delQ0.append(pid[0])
                    E_beam_delQ0.append(pid[1])
                elif i == 1:
                    delE_beam_delQplus1.append(pid[0])
                    E_beam_delQplus1.append(pid[1])
                elif i == 2:
                    delE_beam_delQplus2.append(pid[0])
                    E_beam_delQplus2.append(pid[1])

# ======== Plotting ==========#

    ms = 15 # marker size

    fig0, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 8), sharex=True, sharey=True)
    fig0.suptitle('PID *** {}({},{}){} *** Ebeam = {:.3f} MeV *** dedx = {}'.format(beam, target, ejectile, recoil, beam_energy, dedx_lib),
                  fontname="Times New Roman", fontweight="bold", fontsize = 16)

    ax1.scatter(E_recoil, delE_recoil, label='all recoils', s = 15)
    ax1.scatter(E_beam, delE_beam, c='r', marker='o', s=ms, label=r'leaky beam with E = ({} to {})*E$_b$ '.format(leaky_b_range[0], leaky_b_range[1]))
    ax1.grid()
    ax1.legend(loc='upper left', fontsize=9)

    ax2.scatter(E_recoil_vel, delE_recoil_vel, marker='s', label=r'recoils with v =  $\pm$1.5% v$_{mean}$', s=ms)
    ax2.scatter(E_beam_vel, delE_beam_vel, c='r', marker='s', label='leaky beam with v = $\pm${}% v$_{{mean}}$'.format(100*vel_accept), s = ms)
    ax2.grid()
    ax2.legend(loc='upper left', fontsize=9)

    ax3.scatter(E_recoil_brho, delE_recoil_brho, marker='D', s=ms, label=r'recoils, {}+, B$\rho$ = $\pm$1.5% B$\rho_{{mean}}$'.format(recoil_q))
    ms = 45
    # ax3.scatter(E_beam_delQminus2, delE_beam_delQminus2, c='k', marker='o', s = ms, label=r'leaky beam, {}+, B$\rho$ = $\pm${}% B$\rho_{{mean}}$'.format(beam_q - 2, 100*brho_accept))
    # ax3.scatter(E_beam_delQminus1, delE_beam_delQminus1, c='c', marker='x', s = ms, label=r'leaky beam, {}+, B$\rho$ = $\pm${}% B$\rho_{{mean}}$'.format(beam_q - 1, 100*brho_accept))
    ax3.scatter(E_beam_delQ0, delE_beam_delQ0, c='m', marker='d', s = ms, label=r'leaky beam, {}+, B$\rho$ = $\pm${}% B$\rho_{{mean}}$'.format(beam_q - 0, 100*brho_accept))
    ax3.scatter(E_beam_delQplus1, delE_beam_delQplus1, c='b', marker='+', s = ms, label=r'leaky beam, {}+, B$\rho$ = $\pm${}% B$\rho_{{mean}}$'.format(beam_q + 1, 100*brho_accept))
    ax3.scatter(E_beam_delQplus2, delE_beam_delQplus2, c='g', marker='p', s = ms, label=r'leaky beam, {}+, B$\rho$ = $\pm${}% B$\rho_{{mean}}$'.format(beam_q + 2, 100*brho_accept))
    ax3.set_xlim(E_range[0], E_range[1])
    ax3.set_ylim(delE_range[0], delE_range[1])
    ax3.grid()
    ax3.legend(loc='upper left', fontsize=9)

    fig0.supxlabel('DSSD energy [MeV]')
    fig0.supylabel('IC dE energy [MeV]')

    plt.subplots_adjust(left=0.11, right=0.97, top=0.93, bottom=0.08, hspace=0.0)


    ms = 20 # marker size

    # Converting rad to mrad
    theta_lab = [1e3*x for x in theta_lab]
    theta_lab2 = [1e3*x for x in theta_lab2]
    theta_com = [1e3*x for x in theta_com]

    fig1, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6, 8), sharex=True)
    fig1.suptitle('Kinematics', fontname="Times New Roman", fontweight="bold", fontsize = 16)
    ax1.scatter(theta_lab2, energy_lab, s = ms)
    ax1.set(ylabel = r'$E_{lab}$')
    ax1.grid()
    ax2.scatter(theta_lab, theta_com, s = ms)
    ax2.set(ylabel = r'$\theta_{CM}$ [mrad]')
    ax2.grid()
    ax3.hist(theta_lab, bins=25, histtype='barstacked', label='For isotropic distribution in CoM')
    ax3.set(ylabel = 'counts')
    ax3.grid()
    ax3.legend(loc='upper left', fontsize=10)
    plt.xlabel(r'$\theta_{lab}$ [mrad]', )
    plt.subplots_adjust(left=0.13, right=0.97, top=0.93, bottom=0.08, hspace=0.08)

    plt.show()
