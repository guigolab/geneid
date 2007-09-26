package gphase.model.intransparent.init;

import java.util.Comparator;

public class ASVariationGroup extends ASVariation {

	public ASVariationGroup(ASVariation var) {
		this.trans1= var.trans1;
		this.trans2= var.trans2;
		this.degree= var.degree;
		this.spliceChain1= var.spliceChain1;
		this.spliceChain2= var.spliceChain2;
		this.ssRegionID3UTR= var.ssRegionID3UTR;
		this.ssRegionID5UTR= var.ssRegionID5UTR;
		this.ssRegionIDCDS= var.ssRegionIDCDS;
	}
	
	@Override
	public int hashCode() {
		return getSpliceUniverse()[0].hashCode();
	}
	
	@Override
	public boolean equals(Object obj) {
		Comparator compi= new IdentityComparator();
		return (compi.compare(this, obj)== 0);
	}
}
