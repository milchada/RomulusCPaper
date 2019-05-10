
def sign_angle(v1, v2):

    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    dotp = np.dot(v1_u, v2_u.T)
    return np.sign(dotp.diagonal()) 


        #The colours turned out not to be insightful

        # h1vel = steps[step - (len(t))].halos[0]['Vcom']
        # h2vel = h2ptcls[:1000]['vel'] - h1vel #hm ya how does this relate to central halo before vs after subhalo disruption
            #you can compute this directly with ['vel'] and a dot product between v_rel and pos_rel
        # h2speed = np.linalg.norm(h2vel, axis=1)
        # angle_sign = sign_angle(h2pos, h2vel)
        # h2speed[sign_angle < 0] *= -1
        # norm = colors.Normalize(vmin=200, vmax=800)#h2speed.min(), vmax=h2speed.max())
        # m = cm.ScalarMappable(norm=norm, cmap = cm.viridis)
        # plt.colorbar(m.set_array([]))
        # vel_colors = m.to_rgba(h2speed)
        # print("colours assigned")
