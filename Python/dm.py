#!/usr/bin/env python3


def dm(m , h_eff , i_s , alpha , alpha_l):
    m_is    = np.cross(m,i_s)
    mm_is   = np.cross(m,m_is)
    m_heff  = np.cross(m,h_eff)
    mm_heff = np.cross(m,m_heff)

    dm = - alpha_l * (m_heff + mm_is + alpha * (mm_heff - m_is) )

    return dm
