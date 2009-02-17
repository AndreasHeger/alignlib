//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id$
//--------------------------------------------------------------------------------


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMPL_TOOLKIT_H
#define IMPL_TOOLKIT_H 1

#include "alignlib_fwd.h"
#include "alignlib_interfaces.h"
#include "Toolkit.h"
#include "alignlib.h"

/**
	Basic implementation of toolkit.

   @author Andreas Heger
   @version $Id$: ImplDistor.h,v 1.1.1.1 2002/07/08 21:20:17 heger Exp $
   @short base class of a toolkit
*/


namespace alignlib
{

#define DEFINE_FACTORY(handle,text,make,set,get) \
	public: \
    /** provide a new text object cloned from the default object
     *
     * @return a new @ref text object
     * */ \
    virtual handle make () const { return text->getClone(); } ; \
    /** set the default text object
     */ \
    virtual void set ( const handle & x ) { text = x; }; \
    /** get the default text object
     *
     * @return the default the @ref text object
     */ \
    virtual handle get() const { return text; } ; \
    protected: \
    handle text;

class ImplToolkit: public Toolkit
{
	friend std::ostream & operator<<(std::ostream &output, const Toolkit &);

    // class member functions
 public:
    // constructors and destructors

    /** constructor */
    ImplToolkit() : Toolkit(),
		Alignator( makeAlignatorDPFull( ALIGNMENT_LOCAL, -10, -2) ),
		Fragmentor( alignlib::makeFragmentorDiagonals( Alignator, -10, -2 ) ),
		Alignment( alignlib::makeAlignmentVector() ),
		MultAlignment( alignlib::makeMultAlignment() ),
		MultipleAlignator( alignlib::makeMultipleAlignatorSimple( Alignator ) ),
		Distor( alignlib::makeDistorKimura() ),
		Weightor( alignlib::makeWeightor() ),
		Regularizor( alignlib::makeRegularizor() ),
		LogOddor( alignlib::makeLogOddor()),
		Encoder( alignlib::makeEncoder( Protein20 ) ),
		Treetor( alignlib::makeTreetorDistanceLinkage( Distor ) ),
		Scorer( alignlib::makeScorer() ),
		Iterator2D( alignlib::makeIterator2DFull() ),
		SubstitutionMatrix( alignlib::makeSubstitutionMatrixBlosum62( Encoder ) )
    {};

    /** copy constructor */
    ImplToolkit(const ImplToolkit & src) : Toolkit(src),
    	Alignator(src.Alignator),
    	Fragmentor(src.Fragmentor),
    	Alignment(src.Alignment),
    	MultAlignment(src.MultAlignment),
    	MultipleAlignator(src.MultipleAlignator),
    	Distor( src.Distor),
    	Weightor( src.Weightor),
    	Regularizor( src.Regularizor),
    	LogOddor( src.LogOddor),
    	Encoder( src.Encoder),
    	Treetor( src.Treetor),
    	Scorer( src.Scorer),
    	Iterator2D( src.Iterator2D),
    	SubstitutionMatrix( src.SubstitutionMatrix)
    	{};

    /** destructor */
    virtual ~ImplToolkit () {};

	//------------------------------------------------------------------------------------------------------------
	/** returns a new Toolkit
	 */
	virtual HToolkit getNew() const { return HToolkit(new ImplToolkit) ; };

	/** returns an identical Toolkit
	 */
	virtual HToolkit getClone() const { return HToolkit( new ImplToolkit( *this) ) ; };

    DEFINE_FACTORY( HAlignator, Alignator, makeAlignator, setAlignator, getAlignator);
    DEFINE_FACTORY( HFragmentor, Fragmentor, makeFragmentor, setFragmentor, getFragmentor);
    DEFINE_FACTORY( HAlignment, Alignment, makeAlignment, setAlignment, getAlignment);
    DEFINE_FACTORY( HMultAlignment, MultAlignment, makeMultAlignment, setMultAlignment, getMultAlignment);
    DEFINE_FACTORY( HMultipleAlignator, MultipleAlignator, makeMultipleAlignator, setMultipleAlignator, getMultipleAlignator);
    DEFINE_FACTORY( HDistor, Distor, makeDistor, setDistor, getDistor);
    DEFINE_FACTORY( HWeightor, Weightor, makeWeightor, setWeightor, getWeightor);
    DEFINE_FACTORY( HRegularizor, Regularizor, makeRegularizor, setRegularizor, getRegularizor);
    DEFINE_FACTORY( HLogOddor, LogOddor, makeLogOddor, setLogOddor, getLogOddor);
    DEFINE_FACTORY( HEncoder, Encoder, makeEncoder, setEncoder, getEncoder);
    DEFINE_FACTORY( HTreetor, Treetor, makeTreetor, setTreetor, getTreetor);
    DEFINE_FACTORY( HScorer, Scorer, makeScorer, setScorer, getScorer);
    DEFINE_FACTORY( HIterator2D, Iterator2D, makeIterator2D, setIterator2D, getIterator2D);
    DEFINE_FACTORY( HSubstitutionMatrix, SubstitutionMatrix, makeSubstitutionMatrix, setSubstitutionMatrix, getSubstitutionMatrix);

    virtual void write( std::ostream & output ) const
    {
    	output << "Toolkit - output";
    }
};

HToolkit makeToolkit( const ToolkitType & type )
{
	switch ( type )
	{
	case ProteinAlignment: return HToolkit( new ImplToolkit() ); break;
	case DNAAlignment: return HToolkit( new ImplToolkit() ); break;
	case Genomics: return HToolkit( new ImplToolkit() ); break;
	}

}

}

#endif /* IMPL_TOOLKIT_H */

