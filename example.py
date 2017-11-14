import splat.data
import splat.filters
import splat.gen
import splat.interpol
import splat.scales
from splat import dB2lin as dB

# Create a harmonic scale to use note names rather than frequencies
scale = splat.scales.HarmonicScale()

# Create a triangle wave generator and produce some simple sound with it
triangle = splat.gen.TriangleGenerator()
triangle.filters = [splat.filters.linear_fade]
triangle.run(0.0, 2.8, scale['A-1'], levels=dB(-6.0))

# Create a spline to be used as an envelope (amplitude modulation) in dB
envelope_pts = [(0.0, -100.0), (0.5, -12.0), (1.0, -0.5, 0.0),
                (1.8, -6.0), (2.5, -100.0)]
envelope = splat.interpol.spline(envelope_pts, dB2lin=True)

# Create a sine wave generator and run it with the envelope
sine = splat.gen.SineGenerator()
sine.run(0.0, 2.5, scale['E-1'], levels=envelope.signal)

# Mix the generated fragments together, add some reverb and save the result
master = splat.data.Fragment()
master.mix(triangle.frag, 0.5)
master.mix(sine.frag, 0.75)
splat.filters.reverb(master, splat.filters.reverb_delays())
master.save('example.wav')
