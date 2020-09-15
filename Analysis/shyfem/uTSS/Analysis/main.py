"""
creates gmesh geo file from coastline shapefile
"""
from shapely.geometry import asShape, box as Box, base, mapping, LineString, MultiLineString

import numpy as np
import fiona

def queryCoast(coastLine, bbox):
    print('check0')
    return ((shape, line) for shape, line in
            ((asShape(line["geometry"]), line) for line in coastLine if line["geometry"] is not None)
            if shape.intersects(bbox))


def cutCoast(coastLine, bbox):
    print('check1')
    for shape, line in coastLine:
        shape = shape.intersection(bbox)
        line["geometry"] = mapping(shape.intersection(bbox))
        yield line


def boxCut(file_in, base, box):
    data = fiona.open(base + file_in)
    # data = fiona.open(base + "lines.shp")
    with fiona.open(base + file_in) as coastline:
        # with fiona.open(base + "lines.shp") as coastline:
        meta = coastline.meta
    with fiona.open(base + "output.shp", 'w', **meta) as output:
        for i in range(int(len(data) / 1000) + 1):
            # print i, int(len(data) / 1000) + 1
            for line in cutCoast(queryCoast(data[1000 * i:1000 * (i + 1)], box), box):
                output.write(line)


class Converter:
    def __init__(self):
        self.open = []
        self.closedSizes = [0]
        self.openSizes = [0]
        self.closed = []

    def parse(self, data):
        self.open = []
        self.closed = []
        self.openSizes = [0]
        self.closedSizes = [0]
        for geo in data:
            if geo["type"] == "MultiLineString":
                for c in geo["coordinates"]:
                    self._push(c)
            else:
                self._push(geo["coordinates"])

    def _push(self, coords):
        if coords[0] == coords[-1]:
            self.closed.append(coords[:-1])
            self.closedSizes.append(self.closedSizes[-1] + len(coords) - 1)
        else:
            self.open.append(coords)
            self.openSizes.append(self.openSizes[-1] + len(coords))

    def toGeo(self, data):

        self.parse(data)
        # print type(data)
        coords = np.vstack(self.closed) if len(self.closed) > 0 else np.array([])
        sizes = self.closedSizes
        header = ["// closed curves"]
        indexes = ["IP = newp;", "IL = newl;", "IS = news;", "IF = newf;"]
        points = ["Point(IP + %d) = {%s, %s, 0}; //%d" % (i, c[0], c[1], i)
                  for i, c in enumerate(coords)]
        lines = ["Line(IL + %d) = {IP + %d : IP + %d, IP + %d};" % (
            i, sizes[i], sizes[i + 1] - 1, sizes[i])
                 for i in range(len(sizes) - 1)] + ["//****lastClosed=IL+%d" % (len(sizes) - 1)]

        closedBuffer = header + indexes + points + lines

        coords = np.vstack(self.open) if len(self.open) > 0 else np.array([])
        sizes = self.openSizes
        header = ["// open curves"]

        indexes = ["IP = newp;", "IL = newl;", "IS = news;", "IF = newf;"]
        points = ["Point(IP + %d) = {%s, %s, 0}; //%d" % (i, c[0], c[1], i)
                  for i, c in enumerate(coords)]
        lines = ["Line(IL + %d) = {IP + %d : IP + %d};" % (
            i, sizes[i], sizes[i + 1] - 1)
                 for i in range(len(sizes) - 1)]
        openBuffer = header + indexes + points + lines

        return "\n".join(closedBuffer + openBuffer) + '\n'


class SegmentMerger:
    """
    Merge [Multi]Linestrings when terminations overlap. NO FORK ALOWED
    - closed strings are stored as they are presented
    - isolated strings are stored by terminations
    - if a string as a termination in common with a previously stored string, the two strings are merged together
    """

    def __init__(self):
        self.ends = {}
        self.outs = {}

    @staticmethod
    def join_lines(left, right):
        """ Join segments as 1D numpy array
        :param left: left-side segment
        :param right:  right-side segment
        :return: a new array containing the extended segment
        """
        return np.vstack([left[:-1], right[:]])

    def store(self, item):
        """ Store a segment (by head coordinate)
        :param item: array to store
        """
        self.outs[tuple(item[0])] = item

    def forget(self, item):
        """ Forget a previous stored item
        :param item: the item to forget about
        """
        del self.outs[tuple(item[0])]

    def merge(self, geo):
        """ Merge/Store a new geometry
        :param geo: the segment(s) to merge or store
        :raise: TypeError if geo is a [Multi]LineString
        """
        if geo.type == "MultiLineString":
            [self.merge(i) for i in geo]
        elif geo.type == "LineString":
            coords = np.array(geo.coords)
            start = tuple(coords[0])
            end = tuple(coords[-1])
            if np.all(start == end):
                self.store(coords)
                return
            try:
                other = self.ends[start]
                self.forget(other)
                coords = self.join_lines(other, coords)
            except KeyError:
                pass
            try:
                other = self.ends[end]
                self.forget(other)
                coords = self.join_lines(coords, other)
            except KeyError:
                pass
            self.ends[start] = self.ends[end] = coords
            self.store(coords)

        else:
            raise TypeError("only [Multi]LineString")


def asGeo(shapes, targetBox, simplify=None, updateCb=None):
    updateCb = updateCb if updateCb is not None else lambda *a: None
    features = list(cutCoast(queryCoast(shapes, targetBox), targetBox))
    updateCb("file read")
    shapes = [asShape(i["geometry"]) for i in features[:]]
    updateCb("shaped", len(shapes))
    merger = SegmentMerger()
    for i, shape in enumerate(shapes):
        merger.merge(shape)
        if i % (len(shapes) / 100.) < (i - 1) % (len(shapes) / 100.):
            updateCb("merging", int(i / (len(shapes) / 100.)))
    updateCb("merged")
    shapes = map(LineString, merger.outs.values())
    updateCb("joined")
    shapes = MultiLineString(shapes)
    updateCb("re-joined")
    if simplify is not None:
        shapes = shapes.simplify(simplify, True)
        updateCb("easier now")
    features = [mapping(shape) for shape in shapes]
    conv = Converter()
    return conv.toGeo(features)


def main():

    # path to shapefile
    # mn = fiona.open('/media/milicak/DATA1/datasets/world_bathy/coastlines-WGS84/lines.shp')
    mn = fiona.open('/media/milicak/DATA1/datasets/world_bathy/coastlines-split-4326/lines.shp')

    resol = 0.01

    #base='/Users/scausio/Desktop/'


    # box = (minlon, minlat, maxlon, maxlat)
    # box = Box( -100.0, -40, 60, 80)  # north atlantic
    box = Box( 38, 39, 45, 47)  # Eagean Med
    # box = Box( 22.3, 39, 44, 48)  # Eagean Med

    features = list(cutCoast(queryCoast(mn, box), box))
    # features = list(mn)
    print("file read")
    shapes = [asShape(i["geometry"]) for i in features[:]]

    print("shaped"), len(shapes)
    merger = SegmentMerger()
    for i, shape in enumerate(shapes):
        print(i, shape)
        merger.merge(shape)
        if i % 1000 == 0:
            print (i / 1000)
    print("mehmet")
    shapes = map(LineString, merger.outs.values())
    print("joined")
    shapes = MultiLineString(list(shapes))
    # shapes = MultiLineString(shapes)
    print("re-joined")

    shapes = shapes.simplify(resol, True)
    print("shape simplified")


    # print shapes
    features = [mapping(shape) for shape in shapes]
    conv = Converter()
    prova = conv.toGeo(features)
    print("2geo")
    open("%s.geo"%outName, "w").write(prova)


if __name__ == "__main__":
    base = ""

outName = 'MED_TSS_BS_simplify0_01'
# outName = 'MED_TSS_BS_simplify0_002'
# resol = 0.1
# outName = 'red_sea_simplify%s' % resol

main()

