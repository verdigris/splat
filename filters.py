def linear_fade(frag, duration=0.01):
    fade = min((frag.sample_rate * duration), (len(frag) / 2))
    for i in xrange(int(fade)):
        l = i / fade
        z = ()
        for channel in frag[i]:
            z += ((channel * l),)
        frag[i] = z
        z = ()
        for channel in frag[-i]:
            z += ((channel * l),)
        frag[-i] = z
