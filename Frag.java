import java.util.ArrayList;

/**
 * 
 * @author Nick Fisk
 * Frag is both a storage data type and a node in a graph. This results in a bit of duplicated data
 * but it should be manageable.
 */
public class Frag {
	public String name;
	public String seq;
	public ArrayList<Frag> overlapFrags=new ArrayList<Frag>(); //list of other frags that overlap
	public int amountOverlap;//if this isnt a parent node, this will tell us the number of bases overlapping
	public String finalOverlap;
	public String overlapSeq;//if this isn't a parent node, this will tell us the sequence of the bases overlapping
	
	public Frag(String name, String seq){
		this.name=name;
		this.seq=seq;
		this.finalOverlap="";
	}
}
