/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gribtonetcdfconvertor;

import LRFD.Common.GeoRectangle;
import LRFD.Db.NetCDFFile.NetCDFOperator;
import java.io.File;
import java.io.IOException;
import java.util.*;
import ucar.ma2.DataType;
import ucar.ma2.Index;
import ucar.nc2.*;

/**
 *
 * @author corwin
 */
public class RomsTopLavel 
{
	public class RomsGrid
	{
		NetcdfFile existGrid;
		
		private String lon_rho =  "lon_rho",
				lat_rho = "lat_rho",
				lon_psi = "lon_psi",
				lat_psi = "lat_psi",
				lon_u = "lon_u",
				lat_u = "lat_u",
				lon_v = "lon_v",
				lat_v = "lat_v";
		
		
		
		public class grid
		{
			private double[][] lat;
			private double[][] lon;
			
			public String lon_name, lat_name;
			
			private boolean mCurvilinear = false;
			
			private double[][] loadCoord(String coordinate, boolean isLon) throws IOException
			{
				//NetCDFOperator.Read2DFieldFromFile(inData[t], "Temperature_surface")
				
				Variable var=existGrid.findVariable(coordinate);
				List<Dimension> dims=var.getDimensions();
				
				if (dims.size() != 2)				
					throw new IOException("unexpected lat dimension!");
				
				ucar.ma2.Array varar=var.read();
				
				Index ind=varar.getIndex();
				ind.set(0, 0);
				
				int jn = dims.get(isLon ? 0 : 1 ).getLength();
				int kn = dims.get(isLon ? 1 : 0 ).getLength();
				
				double res[][] = new double[isLon ? jn : kn][isLon ? kn : jn];
				
				for (int j=0; j < jn; j++)
				{
					for (int k=0; k < kn; k++)
					{
						res[isLon ? j : k][isLon ? k : j] =varar.getDouble(ind);
						ind.incr();//set(i,j,k);
						//System.out.print(String.format("%.2f ", SSTar.getDouble(ind)));
					}
				}
			
				return res;
			}
			
			public double []getLat(){ return lat[0];}
			public double []getLon(){ return lon[0];}
			
			public grid (String lonName, String latName) throws IOException
			{
				lon = loadCoord(lon_name = lonName, true);
				lat = loadCoord(lat_name = latName, false);
			}
			
			public  Map<String, Dimension> createDimension(NetcdfFileWriteable cdf) throws IOException
			{
				HashMap<String, Dimension> res = new HashMap();
		
				Dimension lonDim=cdf.addDimension(lon_name, lon[0].length),
						latDim=cdf.addDimension(lat_name, lat[0].length);
		
				Dimension[] lodim=new Dimension[1];
				lodim[0]=lonDim;
				cdf.addVariable(lon_name, DataType.DOUBLE, lodim);
				
				Dimension[] latim=new Dimension[1];
				latim[0] = latDim;
				cdf.addVariable(lat_name, DataType.DOUBLE, latim);
				
				res.put(lat_name, latDim);
				res.put(lon_name, lonDim);
				return res;
			}
		};
		
		public  RomsGrid.grid rho, psi, u, v;
		
		RomsGrid(String road) throws IOException
		{
			existGrid = NetcdfFile.open(road);
			
			rho = new RomsGrid.grid( lon_rho, lat_rho);
			psi = new RomsGrid.grid(lon_psi, lat_psi);
			u = new RomsGrid.grid(lon_u, lat_u);
			v = new RomsGrid.grid(lon_v, lat_v);
		}
		
		public GeoRectangle getRectangle() throws Exception
		{
			return new GeoRectangle( (float)u.getLat()[0] - 4, 
					(float)u.getLat()[u.getLat().length-1] + 4,
					(float)u.getLon()[0] - 4,
					(float)u.getLon()[u.getLon().length-1] + 4);
		}
		
		public Map<String, Dimension> createDimension(NetcdfFileWriteable cdf) throws IOException
		{
			HashMap<String, Dimension> res = new HashMap();
			
			res.putAll(rho.createDimension(cdf));
			res.putAll(psi.createDimension(cdf));
			res.putAll(u.createDimension(cdf));
			res.putAll(v.createDimension(cdf));
			
			return res;
		}
		
		public double[] getParam(String name)
		{
			if (name.equals(lon_rho)) {return rho.getLon();}
			if (name.equals(lat_rho)) {return rho.getLat();}
			if (name.equals(lon_psi)) {return psi.getLon();}
			if (name.equals(lat_psi)) {return psi.getLat();}
			if (name.equals(lon_u)) {return u.getLon();}
			if (name.equals(lat_u)) {return u.getLat();}
			if (name.equals(lon_v)) {return v.getLon();}
			if (name.equals(lat_v)) {return v.getLat();}
			return null;
		}
	};
	
	class RomsVariable
	{
		public String name;
		public String timeName, lonName, latName;
		
		public RomsVariable(String Name, String time, String type)
		{
			name = Name;
			timeName = time;
			lonName = String.format("lon_%s",type);
			latName = String.format("lat_%s",type);
		}
	};
	
	
	public enum VariablesNums 
	{
		Pressure(1),
		Geopotential_height(7),
		Potential_temperature(13),
		Specific_humidity(51),
		Relative_humidity(52),
		SST(11),
		u_wind(33),
		v_wind(34),
		shflux(207),
		dQdSST(208),
		swrad(204),
		svstr(124),
		sustr(125),
		Sensible_heat_flux(122),
		svstr10(127),
		sustr10(128),
		Ground_heat_flux(155),
		Latent_heat_flux(121),
		Ice_thickness(92),
		Ice_concentration_ice1no_ice0(91),
		Land_Surface_Precipitation_Accumulation_LSPA(154),
		Land_cover_land1sea0(81),
		Evaporation(57);	

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
	
	static NetcdfFileWriteable cdf;
	static String dstFile;
	
	RomsGrid grid;
	
	Map<VariablesNums, RomsVariable> resVals;
	
	public RomsTopLavel(String dst_file, String gridFile) throws IOException
	{
		grid = new RomsGrid(gridFile);
		resVals = new HashMap<>(24);
		dstFile = dst_file;
		
		resVals.put(VariablesNums.Pressure, new RomsVariable("Pressure", "time", "u"));
		resVals.put(VariablesNums.Geopotential_height, new RomsVariable("Geopotential_height", "time", "u"));
		
		resVals.put(VariablesNums.Potential_temperature, new RomsVariable("Potential_temperature", "time", "u"));
		resVals.put(VariablesNums.Specific_humidity, new RomsVariable("Specific_humidity", "time", "u"));
		resVals.put(VariablesNums.Relative_humidity, new RomsVariable("Relative_humidity", "time", "u"));
		
		resVals.put(VariablesNums.SST, new RomsVariable("SST", "time", "rho"));//Temperature surfase
		
		resVals.put(VariablesNums.u_wind, new RomsVariable("u_wind", "time", "u"));
		resVals.put(VariablesNums.v_wind, new RomsVariable("v_wind", "time", "v"));
		
		resVals.put(VariablesNums.shflux, new RomsVariable("shflux", "time", "rho")); // sum wave flux (номер указал от балды тк не нашел каой правильный)
		resVals.put(VariablesNums.dQdSST, new RomsVariable("dQdSST", "time", "rho"));
		
		resVals.put(VariablesNums.swrad, new RomsVariable("swrad", "time", "rho")); //Downward_short_wave_flux
		
		resVals.put(VariablesNums.svstr, new RomsVariable("svstr", "time", "v"));//Zonal_momentum_flux
		resVals.put(VariablesNums.sustr, new RomsVariable("sustr", "time", "u")); //Meridional_momentum_flux
		resVals.put(VariablesNums.Sensible_heat_flux, new RomsVariable("Sensible_heat_flux", "time", "u"));
		
		resVals.put(VariablesNums.svstr10, new RomsVariable("svstr10", "time", "v"));
		resVals.put(VariablesNums.sustr10, new RomsVariable("sustr10", "time", "u"));
		
		resVals.put(VariablesNums.Ground_heat_flux, new RomsVariable("Ground_heat_flux", "time", "u"));
		resVals.put(VariablesNums.Latent_heat_flux, new RomsVariable("Latent_heat_flux", "time", "u"));
		
		resVals.put(VariablesNums.Ice_thickness, new RomsVariable("Ice_thickness", "time", "u"));
		resVals.put(VariablesNums.Ice_concentration_ice1no_ice0, new RomsVariable("Ice_concentration_ice1no_ice0", "time", "u"));
		
		resVals.put(VariablesNums.Land_Surface_Precipitation_Accumulation_LSPA, new RomsVariable("Land_Surface_Precipitation_Accumulation_LSPA", "time", "u"));
		resVals.put(VariablesNums.Land_cover_land1sea0, new RomsVariable("Land_cover_land1sea0", "time", "u"));
		resVals.put(VariablesNums.Evaporation, new RomsVariable("Evaporation", "time", "u"));
		createFile();
}
	
	public RomsGrid.grid  getGridForVariable(int varNum)
	{
		switch(varNum)
		{
			case 1:
			case 7:
			case 13:
			case 51:
			case 52:
			case 33:
			case 125:
			case 122:
			case 155:
			case 121:
			case 92:
			case 91:
			case 154:
			case 81:
			case 57:
			case 128:
				return grid.u;
			case 34:
			case 124:
			case 127:
				return grid.v;
			case 11:
			case 204:
			case 207:
			case 208:
				return grid.rho;
			default:
				return null;
		}
	}
	
	public void writeField(Integer fieldNum, double[][][] data)
	{
		NetCDFOperator.writeFieldToNetCDF(dstFile, resVals.get(VariablesNums.getType(fieldNum)).name, data);
	}
	
	public void createFile() throws IOException
	{
		try
		{
			double[] time = {3.0};
			cdf = NetcdfFileWriteable.createNew(dstFile);
			
			Dimension timeDim = cdf.addDimension("time", time.length);
			Map<String,Dimension> Dimensions = grid.createDimension(cdf);

			Dimension[] timeDimArr=new Dimension[1];
			timeDimArr[0]=timeDim;
			cdf.addVariable("time", DataType.DOUBLE, timeDimArr);
			
			Attribute att = new Attribute("missing_value", -9999);
			
			for(Iterator<VariablesNums> i = resVals.keySet().iterator(); i.hasNext();)
			{
				VariablesNums num = i.next();
				RomsVariable var = resVals.get(num);
				Dimension[] todim = new Dimension[3];
				
				todim[0] = timeDim;				
				todim[1] = Dimensions.get(var.latName); 
				todim[2] = Dimensions.get(var.lonName); 
				
				Variable varr = cdf.addVariable(var.name, DataType.DOUBLE, todim);
				
				varr.addAttribute(new Attribute("GRIB_param_number", num.typeValue));
				varr.addAttribute(new Attribute("GRIB_param_name", var.name));
				
				varr.addAttribute(att);
			}
			
			cdf.create();
			cdf.close();
			
			for(Iterator<String> i = Dimensions.keySet().iterator(); i.hasNext();)
			{
				String buf= i.next();
				NetCDFOperator.writeTimeToNetCDF(dstFile, buf, grid.getParam(buf));
			}
			
			NetCDFOperator.writeTimeToNetCDF(dstFile, "time", time);
		}
		catch (Exception er)
		{
			er.printStackTrace();
		}
	}
}
