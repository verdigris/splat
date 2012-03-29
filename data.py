import os
import wave
import _geomusic

# ToDo: make Fragment from file contents
def open_audio_file(file_name, mode):
    ext = file_name.rpartition('.')[2]
    if ext == 'wav':
        return wave.open(file_name, mode)
    raise Exception('Unsupported file format: {0}'.format(ext))


class Fragment(_geomusic.Fragment):
    def resize(self, duration):
        n = int(duration * self.sample_rate)
        self._resize(n)
        return n

    def save_to_file(self, file_name, sample_width):
        print("sample width: {0}, channels: {1}, length: {2}".format(
                sample_width, self.channels, len(self)))
        f = wave.open(file_name, 'w')
        f.setnchannels(self.channels)
        f.setsampwidth(sample_width)
        f.setframerate(self.sample_rate)
        f.setnframes(len(self))
        f.writeframes(self.as_bytes(sample_width))
        f.close()
