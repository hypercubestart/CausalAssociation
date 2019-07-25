import numpy

class SingleCell:
    def __init__(self, variant_count, regulon_count):
        self.variants = numpy.zeros(variant_count)
        self.regulon_activity = numpy.zeros(regulon_count)

    def set_variant(self, index, value):
        self.variants[index] = value

    def set_regulon(self, index, value):
        self.regulon_activity[index] = value

    def get_variant(self, index):
        return self.variants[index]

    def get_regulon(self, index):
        return self.regulon_activity[index]

    def copy_regulon(self):
        return self.regulon_activity.copy()
