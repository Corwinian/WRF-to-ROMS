/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gribtonetcdfconvertor;

/**
 *
 * @author corwin
 */
public enum VariablesNums 
{
	Pressure(1),
	Geopotential_height(7),
	Potential_temperature(13),
	Specific_humidity(51),
	Relative_humidity(52),
	SST(11),
	Temperature(12),
	u_wind(33),
	v_wind(34),
	swflux(309),
	shflux(207),
	dQdSST(208),
	swrad (206),
	dsrad (204),//Incoming surface shortwave radiation — time-averaged
	dlrad(205),//Incoming surface longwave radiation - time-averaged
	usrad(211),//Outgoing surface shortwave radiation - time-averaged	
	ulrad(212),//Outgoing surface longwave radiation – time-averaged
	svstr(124),
	sustr(125),
	Sensible_heat_flux(122),	
	Ground_heat_flux(155),
	Latent_heat_flux(121),
	Ice_thickness(92),
	Ice_concentration_ice1no_ice0(91),
	Land_Surface_Precipitation_Accumulation_LSPA(154),
	Land_cover_land1sea0(81),
	Evaporation(57),
	Surface_moment_flux(172),//Surface momentum flux — time-averaged
	SSS(175)//
	;	

	private Integer typeValue;
	private VariablesNums(Integer type) {typeValue = type;}
	public Integer getTypeValue() {return typeValue;}

	static public VariablesNums  getType(Integer pType) 
	{
		for (VariablesNums  type: VariablesNums.values()) 
		{
			if (type.getTypeValue().equals(pType)) 
				return type;
		}
		throw new RuntimeException("unknown type");
	}
}

