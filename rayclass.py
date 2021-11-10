from sage.groups.abelian_gps.values import AbelianGroupWithValues_class, AbelianGroupWithValuesElement

class RayIdealClass(AbelianGroupWithValuesElement):
    """NB: this might be merged into FractionalIdealClass at some
    point, but for now that doesn't support ideals of orders.

    """

    def __init__(self, parent, element,):
        if element is None:
            element = parent._ideal_log(ideal)
        AbelianGroupWithValuesElement.__init__(self, parent, element, ideal)

    def _mul_(self, other):
        m = AbelianGroupElement._mul_(self, other)
        m._value = 1            # todo
        return m

    def _div_(self, other):
        d = AbelianGroupElement._mul_(self, other)
        m._value = 1             # todo
        return m

    def __pow__(self, n):
        # We use MonoidElement's __pow__ routine, since that does
        # repeated squaring, and hence the ideal gets reduced as
        # we go along; actually computing self._value ** n would
        # be disastrous.
        n = n % self.order()
        return MonoidElement.__pow__(self, n)

    def inverse(self):
        r"""
        Return the multiplicative inverse of this ideal class.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 - 3*x + 8); G = K.class_group()
            sage: G(2, a).inverse()
            Fractional ideal class (2, a^2 + 2*a - 1)
            sage: ~G(2, a)
            Fractional ideal class (2, a^2 + 2*a - 1)
        """
        m = AbelianGroupElement.inverse(self)
        m._value = 1             # todo
        return m
    __invert__ = inverse

    def _repr_(self):
        if self.is_principal():
            return 'Trivial principal ray ideal class'
        return 'Fractional ray ideal class %s'%self._value._repr_short()

    def is_principal(self):
        """
        Test if an ray ideal class is principal
        """
        return(1)


class RayClassGroup():
    """
    Class of ray class groups
    """
    # Element = RayIdealClass
    def __init__(self, number_field, modulus):
        r"""Create ray class group $\mathrm{Cl}_m(F)$ where $F$ is a
        number field and $m$ is the modulus.
        
        Input:
        - Number Field F
        - Modulus m, which is either an ideal, or a list consisting of 
          - an integral ideal I
          - a list of length r1 corresponding to the real infinite places of F

        """

        self.number_field = number_field
        self.modulus = modulus
        # AbelianGroupWithValues_class.__init__(self, gens_orders, names, gens,
        #                                       values_group=number_field.ideal_monoid())
        # self.modulus = modulus
        # self.order = 

    # def _element_constructor_(self, *args, **kwds):
    #     I = self._number_field.ideal(*args, **kwds)
    #     return self.element_class(self, None, I)

    # def __iter__(self):        
    def group(self):

        """
        Gives the underlying additive group
        """
        F = self.number_field
        m = self.modulus
        P = F.pari_polynomial()

        # r1 = Integer(F.signature()[0])
        # Initialise pari field
        bnf = pari("bnfinit({})".format(P))
        # For the time being, let's assume m is principal
        # decompose m in pari
        mPari = pari("idealhnf({},{})".format(bnf, m))

        # compute Cl_m^+
        L = pari("bnrinit({},[{},{}],1).cyc".format(bnf, mPari, [1]*F.signature()[0]))
        if len(L) == 0:
            return(AdditiveAbelianGroup([1]))
        return(AdditiveAbelianGroup(L))


# class RayClassCharacterGroup(AbelianGroupWithValues_class, UniqueRepresentation):
#     """
#     """
#     Element = RayClassCharacter
#     def __init__(self, number_field, discriminant=1):
#         """
#         Create character group
#         """
#         self.number_field = number_field
#         self.ray_class_group =
#         self.modulus = 

# class RayClassCharacter (AbelianGroupWithValuesElement):
#     """
#     """

#     def __init__(self, parent,element,):
#         if element is None:
#             element = parent._ideal_log(ideal)
#         AbelianGroupWithValuesElement.__init__(self, parent, element, ideal)

#     def _mul_(self, other):
#         m = AbelianGroupElement._mul_(self, other)
#         m._value =              # todo
#         return m

#     def _div_(self, other):
#         d = AbelianGroupElement._mul_(self, other)
#         m._value =              #
#         return m

#     def __pow__(self, n):
#         # We use MonoidElement's __pow__ routine, since that does
#         # repeated squaring, and hence the ideal gets reduced as
#         # we go along; actually computing self._value ** n would
#         # be disastrous.
#         n = n % self.order()
#         return MonoidElement.__pow__(self, n)

#     def inverse(self):
#         r"""
#         Return the multiplicative inverse of this ideal class.

#         EXAMPLES::

#             sage: K.<a> = NumberField(x^3 - 3*x + 8); G = K.class_group()
#             sage: G(2, a).inverse()
#             Fractional ideal class (2, a^2 + 2*a - 1)
#             sage: ~G(2, a)
#             Fractional ideal class (2, a^2 + 2*a - 1)
#         """
#         m = AbelianGroupElement.inverse(self)
#         m._value =              # 
#         return m
#     __invert__ = inverse
