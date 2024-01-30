import sys

## classes and helper functions
class TaxonLevel:
    ## make the class behave
    def __init__(self, tax_lvl):

        self.lvl = tax_lvl
        if len(tax_lvl) == 2:
            self.val = int(tax_lvl[1])
        else:
            self.val = 0

    ## dress up the class properly
    def __repr__(self):
        return f"TaxonLevel({self.lvl})"
    def __str__(self):
        return self.lvl

    ## compare TaxonLevels
    def __eq__(self, other):
        if not isinstance(other, TaxonLevel):
            raise TypeError(f"{other} not of type TaxonLevel.")
        return self.val == other.val

    def __gt__(self, other):
        if not isinstance(other, TaxonLevel):
            raise TypeError(f"{other} not of type TaxonLevel.")
        return self.val > other.val

    def __ls__(self, other):
        if not isinstance(other, TaxonLevel):
            raise TypeError(f"{other} not of type TaxonLevel.")
        return self.val < other.val

    def __ge__(self, other):
        if not isinstance(other, TaxonLevel):
            raise TypeError(f"{other} not of type TaxonLevel.")
        return self.val >= other.val

    def __le__(self, other):
        if not isinstance(other, TaxonLevel):
            raise TypeError(f"{other} not of type TaxonLevel.")
        return self.val <= other.val

    ## do math on TaxonLevels
    def __add__(self, value):
        if isinstance(value, int) or isinstance(value, float):
            return TaxonLevel("S"+str(self.val + int(value)))
        if isinstance(value, TaxonLevel):
            updated_val = self.val + value.val
            return TaxonLevel("S"+str(updated_val))

    def __sub__(self, value):
        if isinstance(value, int) or isinstance(value, float):
            work_factor = self.val - value
            if work_factor < 0:
                sys.stderr.write("Too far, returning to S")
                return TaxonLevel("S")
            if work_factor == 0:
                return TaxonLevel("S")
            else:
                return TaxonLevel("S"+str(work_factor))

        if isinstance(value, TaxonLevel):
            work_factor = self.val - value.val
            if work_factor < 0:
                sys.stderr.write("Too far, returning to S")
            else:
                return work_factor

    ## find logical paths between TaxonLevels
    def to(self, other):
        if not isinstance(other, TaxonLevel):
            raise TypeError(f"{other} not of type TaxonLevel.")

        return ["S"+str(i) for i in range(self.val+1, other.val)]

## find indices of a target element in code
def find_locs(target, input_list):
    return [i for i, n in enumerate(input_list) if n == target]

## find closest parent
def find_parent(input_list, parent):
    for i, e in reversed(list(enumerate(input_list))):
        if e == parent:
            return i