#!/usr/bin/env python3


def rk4(m , h_eff , i_s , dt):
    k = np.zeros((4,3))
    mm = np.zeros((3,3))
    m_new = np.zeros((1,3))

    # Step 0
    k[0 , :] = dm(m , h_eff , i_s)
    mm[0 , :] = m + k[0 , :] * dt / 2

    # Step 1
    k[1 , :] = dm(mm[0 , :] , h_eff , i_s)
    mm[1 , :] = m + k[1 , :] * dt / 2

    # Step 2
    k[2 , :] = dm(mm[1 , :] , h_eff , i_s)
    mm[2 , :] = m + k[2 , :] * dt

    # Step 3
    k[3 , :] = dm(mm[2 , :] , h_eff , i_s)

    # Return new 'm'
    m_new = (k[0 , :] + 2 * (k[1 , :] + k[2 , :]) + k[3 , :])* dt / 6
    return m_new
