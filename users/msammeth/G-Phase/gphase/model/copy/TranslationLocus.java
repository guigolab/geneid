package gphase.model.copy;

public class TranslationLocus extends Gene {

	static final long serialVersionUID = 8933737248601221991L;
	
	Gene gene= null;
	
	public TranslationLocus(Gene ge, Transcript[] trpts) {
		super(ge.getSpecies(), ge.getID());
		this.gene= ge;
		this.transcripts= trpts;
	}
	
	
}
