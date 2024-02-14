import math
import random

# Functions kinematics, cm2lab, and cm2lab_theta are translated from fortran code STOPIT.

def kinematics(m1:float, m2:float, m3:float, m4:float, exc3:float, eb:float, t_lab:float):
    m3 = m3 + exc3

    w2l = eb + m2
    g2l = w2l / m2
    b2l = (1.0 - 1.0 / (g2l * g2l)) ** 0.5
    p2l = g2l * b2l * m2

    w1l = m1
    wtl = w1l + w2l

    b = p2l / (w2l + m1)
    g = 1.0 / (1.0 - b * b) ** 0.5

    wtc = wtl / g

    e3l = 0

    if wtc < (m3 + m4):
        return e3l

    elif wtc > (m3 + m4):
        solution = []
        w3c = (wtc * wtc + m3 * m3 - m4 * m4) / (2 * wtc)
        p3c = (w3c * w3c - m3 * m3) ** 0.5
        g3c = w3c / m3
        b3c = (1.0 - 1.0 / (g3c * g3c)) ** 0.5
        k3c = b / b3c

        w4c = (wtc * wtc + m4 * m4 - m3 * m3) / (2 * wtc)
        p4c = (w4c * w4c - m4 * m4) ** 0.5
        if m4 == 0:
            b4c = 1
        else:
            g4c = w4c / m4
            b4c = (1.0 - 1.0 / (g4c * g4c)) ** 0.5

        k4c = b / b4c

        t3l = t_lab
        t3c = 0
        xx = (m1 - m4) * (m1 + m4) + (m2 - m3) * (2 * m1 + m2 - m3) + 2 * eb * (m1 - m3)
        aa = math.cos(t3l) * math.cos(t3l) + g * g * math.sin(t3l) * math.sin(t3l)
        bb = (g * math.sin(t3l)) ** 2 * 2.0 * k3c
        cc = 4. * math.cos(t3l) * math.cos(t3l) * (math.cos(t3l) * math.cos(t3l) + g * g * math.sin(t3l)
                                                   * math.sin(t3l) * xx * (g3c / g + 1.0) / (
                                                           (g3c * g3c - 1.0) * (2.0 * m3 * wtl)))

        solution = []
        tmax = 0
        if xx > 0:  # theta max doesn't exist and there is only 1 solution.
            if math.cos(t3l) >= 0.0:
                t3c = math.acos((-bb + (cc) ** 0.5) / (2.0 * aa))
            elif math.cos(t3l) < 0:
                t3c = math.acos((-bb - (cc) ** 0.5) / (2.0 * aa))
            solution.append(cm2lab(m3, t3c, w3c, b3c, m4, w4c, b4c, k4c, g, b)[0])

        elif xx <= 0:  # theta max exists and there are 2 solutions.
            xa = -g * xx * (g3c + g) / (2. * m3 * wtl * (g3c * g3c - 1.0))
            tmax = math.acos((xa / (xa + 1.0)) ** 0.5)
            if t3l <= tmax:
                t3c = math.acos((-bb - (cc) ** 0.5) / (2.0 * aa))
                solution.append(cm2lab(m3, t3c, w3c, b3c, m4, w4c, b4c, k4c, g, b)[0])
                t3c = math.acos((-bb + (cc) ** 0.5) / (2.0 * aa))
                solution.append(cm2lab(m3, t3c, w3c, b3c, m4, w4c, b4c, k4c, g, b)[0])

        solution.append(tmax)

    return solution


def cm2lab(m3, t3c, w3c, b3c, m4, w4c, b4c, k4c, g, b):
    t4c = -math.pi + t3c
    t4l = math.atan(math.sin(t4c) / (g * (math.cos(t4c) + k4c)))
    if t4c * t4l < 0.0:
        t4l = t4l - math.pi
    w3l = w3c * g * (1.0 + b * b3c * math.cos(t3c))
    e3l = w3l - m3
    w4l = w4c * g * (1.0 + b * b4c * math.cos(t4c))
    e4l = w4l - m4
    solution = [e3l, t4l, e4l]
    return solution


def cm2lab_theta(t3c, m1, m2, m3, m4, exc3, eb):
    m3 = m3 + exc3
    w2l = eb + m2
    g2l = w2l / m2
    b2l = (1.0 - 1.0 / (g2l * g2l)) ** 0.5
    p2l = g2l * b2l * m2
    w1l = m1
    wtl = w1l + w2l
    b = p2l / (w2l + m1)
    g = 1.0 / (1.0 - b * b) ** 0.5
    wtc = wtl / g
    e3l = 0
    if wtc < (m3 + m4):
        return e3l
    elif wtc > (m3 + m4):
        w3c = (wtc * wtc + m3 * m3 - m4 * m4) / (2 * wtc)
        g3c = w3c / m3
        b3c = (1.0 - 1.0 / (g3c * g3c)) ** 0.5
    w3l = w3c * g * (1.0 + b * b3c * math.cos(t3c))
    t3l = math.asin( math.sin(t3c) * ((w3c**2-m3**2)/(w3l**2 - m3**2))**0.5)
    return t3l


def generate_random_theta_com():
    rs = 0
    while rs > 1 or rs < 0.01:
        x = random.uniform(-1,1)
        y = random.uniform(-1,1)
        z = random.uniform(-1,1)
        rs = (x**2 + y**2 + z**2)**0.5
    el = (abs(x**2 + y**2))**0.5
    if abs(el) < 1e-9: el += 3e-9
    if abs(z) < 1e-9: z += 3e-9
    th = math.acos(z/rs)
    return th


def get_brho(energy, mass, charge):
    total_energy = energy*1e-3 + mass
    pc = (total_energy**2 - mass**2)**0.5
    return 3.33564*pc/charge


def get_vel(energy, mass):
    c = 29.9792 # cm/ns
    total_energy = energy * 1e-3 + mass
    return c*(1-(mass/total_energy)**2)**0.5

