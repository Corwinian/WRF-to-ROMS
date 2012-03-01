/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gribtonetcdfconvertor;

import LRFD.Common.GeoRectangle;
import LRFD.Db.NetCDFFile.NetCDFOperator;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import ucar.ma2.DataType;
import ucar.ma2.Index;
import ucar.nc2.*;

/**
 *
 * @author corwin
 */
public final class RomsTopLavel 
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
				lat_v = "lat_v",

				srt_time = "srt_time",
				sst_time = "sst_time",
				sss_time = "sss_time",
				sms_time = "sms_time",
				shf_time = "shf_time",
				swf_time = "swf_time",
				g_time = "g_time";
		
		
		
		public class grid
		{
			private double[][] lat;
			private double[][] lon;
			
			public String lon_name, lat_name;
			
			private boolean mCurvilinear = false;
			
			private double[][] loadCoord(String coordinate, boolean isLon) throws IOException
			{
				//NetCDFOperator.Read2DFieldFromFile(inData[t], "Temperature_surface")
				
				Variable var = existGrid.findVariable(coordinate);
				List<Dimension> dims = var.getDimensions();
				
				if (dims.size() != 2)				
					throw new IOException("unexpected lat dimension!");
				
				ucar.ma2.Array varar = var.read();
				
				Index ind = varar.getIndex();
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
		
				Dimension lonDim = cdf.addDimension(lon_name, lon[0].length),
						latDim = cdf.addDimension(lat_name, lat[0].length);
		
				Dimension[] lodim = new Dimension[1];
				lodim[0] = lonDim;
				cdf.addVariable(lon_name, DataType.DOUBLE, lodim);
				
				Dimension[] latim = new Dimension[1];
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
		
		public Map<String, Dimension> createDimension(NetcdfFileWriteable cdf, Integer timeLen) throws IOException
		{
			HashMap<String, Dimension> res = new HashMap();
			
			res.putAll(rho.createDimension(cdf));
			res.putAll(psi.createDimension(cdf));
			res.putAll(u.createDimension(cdf));
			res.putAll(v.createDimension(cdf));
			
			String[] times = {srt_time, sst_time, sss_time, sms_time, shf_time, swf_time, g_time};
		
			for (int i=0; i < times.length; ++i)
			{
				Dimension timeDim = cdf.addDimension(times[i], timeLen);
			
				Dimension[] tidim = new Dimension[1];
				tidim[0] = timeDim;
				cdf.addVariable(times[i], DataType.DOUBLE, tidim);
				
				res.put(times[i], timeDim);
			}
			
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
			timeName = String.format("%s_time", time);
			lonName = String.format("lon_%s", type);
			latName = String.format("lat_%s", type);
		}
	};
	
	
	static NetcdfFileWriteable cdf;
	static String dstFile;
	
	RomsGrid grid;
	
	Map<VariablesNums, RomsVariable> resVals;
	
	public Map<VariablesNums, RomsVariable> getResVals() {return resVals;}
	
	public RomsTopLavel(String dst_file, String gridFile, double  []time) throws IOException
	{
		grid = new RomsGrid(gridFile);
		resVals = new HashMap<>(24);
		dstFile = dst_file;
		
		//TODO:
		resVals.put(VariablesNums.Pressure, new RomsVariable("Pressure", "g", "u"));
		resVals.put(VariablesNums.Geopotential_height, new RomsVariable("Geopotential_height", "g", "u"));
		
		resVals.put(VariablesNums.Potential_temperature, new RomsVariable("Potential_temperature", "g", "u"));
		resVals.put(VariablesNums.Specific_humidity, new RomsVariable("Specific_humidity", "g", "u"));
		resVals.put(VariablesNums.Relative_humidity, new RomsVariable("Relative_humidity", "g", "u"));
		
		resVals.put(VariablesNums.SST, new RomsVariable("SST", "sst", "rho"));//Temperature surfase
		
		resVals.put(VariablesNums.u_wind, new RomsVariable("u_wind", "g", "u"));
		resVals.put(VariablesNums.v_wind, new RomsVariable("v_wind", "g", "v"));
		
		resVals.put(VariablesNums.shflux, new RomsVariable("shflux", "shf", "rho")); // sum wave flux (номер указал от балды тк не нашел каой правильный)
		resVals.put(VariablesNums.dQdSST, new RomsVariable("dQdSST", "sst", "rho"));
		
		resVals.put(VariablesNums.swrad, new RomsVariable("swrad", "srt", "rho")); //Downward_short_wave_flux
		
		resVals.put(VariablesNums.svstr, new RomsVariable("svstr", "sms", "v"));//Zonal_momentum_flux
		resVals.put(VariablesNums.sustr, new RomsVariable("sustr", "sms", "u")); //Meridional_momentum_flux
		//FEXME: подумать над названием 
		resVals.put(VariablesNums.Sensible_heat_flux, new RomsVariable("Sensible_heat_flux", "g", "u"));
		
		resVals.put(VariablesNums.Ground_heat_flux, new RomsVariable("Ground_heat_flux", "g", "u"));
		resVals.put(VariablesNums.Latent_heat_flux, new RomsVariable("Latent_heat_flux", "g", "u"));
		
		resVals.put(VariablesNums.Ice_thickness, new RomsVariable("Ice_thickness", "g", "u"));
		resVals.put(VariablesNums.Ice_concentration_ice1no_ice0, new RomsVariable("Ice_concentration_ice1no_ice0", "g", "u"));
		
		resVals.put(VariablesNums.Land_Surface_Precipitation_Accumulation_LSPA, new RomsVariable("Land_Surface_Precipitation_Accumulation_LSPA", "g", "u"));
		resVals.put(VariablesNums.Land_cover_land1sea0, new RomsVariable("Land_cover_land1sea0", "g", "u"));
		resVals.put(VariablesNums.Evaporation, new RomsVariable("Evaporation", "g", "u"));
		createFile(time);
	}
	
	public RomsGrid.grid  getGridForVariable(VariablesNums varNum)
	{
		switch(varNum)
		{
			case Pressure:
			case Geopotential_height:
			case Potential_temperature:
			case Specific_humidity:
			case Relative_humidity:
			case u_wind:
			case sustr:
			case Sensible_heat_flux:
			case Ground_heat_flux:
			case Latent_heat_flux:
			case Ice_thickness:
			case Ice_concentration_ice1no_ice0:
			case Land_Surface_Precipitation_Accumulation_LSPA:
			case Land_cover_land1sea0:
			case Evaporation:
				return grid.u;
			case v_wind:
			case svstr:
				return grid.v;
			case SST:
			case swrad:
			case shflux:
			case dQdSST:
				return grid.rho;
			default:
				return null;
		}
	}
	
	public void writeField(VariablesNums fieldNum, double[][][] data)
	{
		NetCDFOperator.writeFieldToNetCDF(dstFile, resVals.get(fieldNum).name, data);
	}
	
	public void createFile(double []time) throws IOException
	{
		try
		{
			cdf = NetcdfFileWriteable.createNew(dstFile);
			
			Map<String,Dimension> Dimensions = grid.createDimension(cdf, time.length);

			Attribute att = new Attribute("missing_value", -9999);
			
			for(Iterator<VariablesNums> i = resVals.keySet().iterator(); i.hasNext();)
			{
				VariablesNums num = i.next();
				RomsVariable var = resVals.get(num);
				Dimension[] todim = new Dimension[3];
				
				todim[0] = Dimensions.get(var.timeName);				
				todim[1] = Dimensions.get(var.latName); 
				todim[2] = Dimensions.get(var.lonName); 
				
				Variable varr = cdf.addVariable(var.name, DataType.DOUBLE, todim);
				
	 			varr.addAttribute(new Attribute("GRIB_param_number", num.getTypeValue()));
				varr.addAttribute(new Attribute("GRIB_param_name", var.name));
				
				varr.addAttribute(att);
			}
			
			cdf.create();
			cdf.close();
			
			for(Iterator<String> i = Dimensions.keySet().iterator(); i.hasNext();)
			{
				String buf = i.next();
				NetCDFOperator.writeTimeToNetCDF(dstFile, buf,
						buf.substring(4).equals("time") ||  buf.substring(2).equals("time") 
						?  time : grid.getParam(buf));
			}
		}
		catch (Exception er)
		{
			er.printStackTrace();
		}
	}
}
