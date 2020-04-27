import bisect
import os
import os.path

DIRNAME = '/media/vineetb/t5-vineetb/dsprint/in/hgdownload/hg19/phastCons100way/hg19.100way.phastCons'
FILENAME = 'chr4.phastCons100way.wigFix'


class WigFix:

    def __init__(self, filename):
        self.filename = filename
        positions = []
        last_position = 0
        values = {}

        with open(filename, 'r') as f:
            line = f.readline()
            while line != '':
                if line.startswith('fixedStep'):
                    last_position = int(line.split(' ')[2].split('=')[1])
                    values[last_position] = []
                    positions.append(last_position)
                else:
                    values[last_position].append(float(line))
                line = f.readline()

        self.positions = positions
        self.values = values

    def __getitem__(self, item):
        position = bisect.bisect_right(self.positions, item)
        # The return position from bisect_right is the insert position
        # This is 0 for elements < the first that we have, 1 between [<first>, <second>)
        # Subtract 1 to get the index where we can start our forward search
        if position < 1:
            return None
        start_position = self.positions[position - 1]
        try:
            return self.values[start_position][item - start_position]
        except IndexError:
            return None


if __name__ == '__main__':

    o = WigFix(os.path.join(DIRNAME, FILENAME))
    print(o[10528])
    print(o[11528])
    print(o[191041279])
    print(o[191041280])
    print(o[191041281])
    print(o[192041281])
