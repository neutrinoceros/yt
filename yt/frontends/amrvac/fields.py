"""
AMRVAC-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer

code_density = "code_density"
code_velocity = "code_length / code_time"
code_momentum = "code_density * code_velocity / code_time"
code_temperature = "code_temperature"
code_pressure = "code_density * code_temperature"
code_energy = "code_specific_energy"
code_magnetic = "code_magnetic"

# adiabatic constant for non HD/MHD datasets, used in the EoS for pressure
adiab_cte = 1.0

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class AMRVACFieldInfo(FieldInfoContainer):
    known_other_fields = (
        # Everything in AMRVAC is normalised to dimensionless units, so set units to ""
        ("rho", (code_density, ["density"], r"$\rho$")),
        ("m1", (code_momentum, ["momentum_1"], r"$m_1$")),
        ("m2", (code_momentum, ["momentum_2"], r"$m_2$")),
        ("m3", (code_momentum, ["momentum_3"], r"$m_3$")),
        ("e", (code_energy, ["energy"], r"$e$")),
        ("b1", (code_magnetic, [], r"$B_1$")),
        ("b2", (code_magnetic, [], r"$B_2$")),
        ("b3", (code_magnetic, [], r"$B_3$"))
    )

    known_particle_fields = ()

    def __init__(self, ds, field_list):
        super(AMRVACFieldInfo, self).__init__(ds, field_list)

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases

        # add primitive variables as custom fields
        def _get_ekin(data):
            ekin = 0.5 * data["amrvac", "m1"]**2 / data["amrvac", "rho"]
            if self.ds.dimensionality > 1:
                ekin = ekin + 0.5 * data["amrvac", "m2"] ** 2 / data["amrvac", "rho"]
            if self.ds.dimensionality > 2:
                ekin = ekin + 0.5 * data["amrvac", "m3"] ** 2 / data["amrvac", "rho"]
            ekin.units = code_energy
            return ekin

        def _get_emag(data):
            emag = 0.5 * data["amrvac", "b1"]**2
            if self.ds.dimensionality > 1:
                emag = emag + 0.5 * data["amrvac", "b2"]**2
            if self.ds.dimensionality > 2:
                emag = emag + 0.5 * data["amrvac", "b3"]**2
            emag.units = code_energy
            return emag

        def _pressure(field, data):
            # (M)HD datasets
            if ("amrvac", "e") in self.field_list:
                # MHD dataset
                if ("amrvac", "b1") in self.field_list:
                    pres = (self.ds.gamma - 1) * (data["amrvac", "e"] - _get_ekin(data) - _get_emag(data))
                # Non-MHD dataset
                else:
                    pres = (self.ds.gamma - 1) * (data["amrvac", "e"] - _get_ekin(data))
            else:
                # TODO: define unit_pressure in case there is no energy equation (right now it's code_pressure above)
                # Non (M)HD datasets, in this case an EoS is used for the pressure
                pres = adiab_cte * data["amrvac", "rho"] ** self.ds.parameters.gamma
            pres.units = code_pressure
            return pres

        def _temperature(field, data):
            return data["amrvac", "pressure"] / data["amrvac", "rho"]

        def _velocity1(field, data):
            return data["amrvac", "m1"] / data["amrvac", "rho"]

        def _velocity2(field, data):
            return data["amrvac", "m2"] / data["amrvac", "rho"]

        def _velocity3(field, data):
            return data["amrvac", "m3"] / data["amrvac", "rho"]

        def _total_energy(field, data):
            etot = _get_ekin(data)
            if ("amrvac", "b1") in self.field_list:
                etot = etot + _get_emag(data)
            return etot

        # pressure field, add this first to calculate temperature after
        self.add_field(("amrvac", "pressure"), sampling_type="cell",
                       function=_pressure, units='')
        self.alias(("gas", "pressure"), ("amrvac", "pressure"), units='')

        # temperature field
        self.add_field(("amrvac", "temperature"), sampling_type="cell",
                       function=_temperature, units='')
        self.alias(("gas", "temperature"), ("amrvac", "temperature"), units='')

        # velocity fields
        self.add_field(("amrvac", "velocity_1"), sampling_type="cell",
                       function=_velocity1, units="")
        self.alias(("amrvac", "v1"), ("amrvac", "velocity_1"), units="")
        if ("amrvac", "m2") in self.field_list:
            self.add_field(("amrvac", "velocity_2"), sampling_type="cell",
                           function=_velocity2, units="")
            self.alias(("amrvac", "v2"), ("amrvac", "velocity_2"), units="")
        if ("amrvac", "m3") in self.field_list:
            self.add_field(("amrvac", "velocity_3"), sampling_type="cell",
                           function=_velocity3, units="")
            self.alias(("amrvac", "v3"), ("amrvac", "velocity_3"), units="")

        # total energy
        self.add_field(("amrvac", "total_energy"), sampling_type="cell",
                       function=_total_energy, units="")
        self.alias(("gas", "total_energy"), ("amrvac", "total_energy"), units="")

        setup_magnetic_field_aliases(self, "amrvac", ["mag%s" % ax for ax in "xyz"])


    def setup_particle_fields(self, ptype):
        super(AMRVACFieldInfo, self).setup_particle_fields(ptype)
