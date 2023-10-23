import java.util.ArrayList;

/**
 * Contig is a simple container class to hold sequences that
 * have been merged together into a contig. 
 * @author Nick Fisk
 *
 */

public class Contig{
	//list of all the component fragments
	public ArrayList<String> frags=new ArrayList<String>();
	//consensus sequence
	public String consensus;
	
	/*
	 * constructor for a contig instance.
	 */
	public Contig(){
		this.consensus="";
	
	}
}
