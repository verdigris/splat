def linear_fade(frag, duration=0.01):
    fade = min((frag.sample_rate * duration), (len(frag) / 2))
    for i in range(int(fade)):
        l = i / fade
        for j in (i, -i):
            frag[j] = tuple((s * l) for s in frag[j])
