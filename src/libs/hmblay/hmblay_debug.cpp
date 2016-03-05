#include "hmblay_debug.hpp"

#ifndef NDEBUG
using namespace HMBlay;

void Debug::info_extpath(const Impl::ExtPath& pth){
	std::cout<<"+++ Extended path at "<<&pth<<". "<<pth.size()<<" edges. ";
	if (pth.is_closed()) std::cout<<"Closed. "<<"Area = "<<HMCont2D::Contour::Area(pth)<<std::endl;
	else std::cout<<"Open."<<std::endl;
	std::cout<<"+++ Length = "<<pth.length()<<std::endl;
	auto op = pth.ordered_points();
	for (int i=0; i<op.size(); ++i){
		printf("Point %p:  (%10.6f, %10.6f)", op[i], op[i]->x, op[i]->y);
		if (i<pth.size())
			printf("   --> Edge %p", pth.edge(i));
		else
			printf("   --> Edge %p", pth.edge(i-pth.size()));
		
		double an = pth.ext_data[i].angle/M_PI*180;
		printf("  angle = %6.2f;", an);

		double nangle = ToAngle(atan2(pth.ext_data[i].normal.y, pth.ext_data[i].normal.x));
		nangle = nangle/M_PI*180;
		printf("  normal angle = %6.2f;", nangle);
		
		printf("   ");
		switch (pth.ext_data[i].tp){
			case Impl::CornerTp::NO: printf("No"); break;
			case Impl::CornerTp::ZERO: printf("Zero"); break;
			case Impl::CornerTp::STRAIGHT: printf("Straight"); break;
			case Impl::CornerTp::RIGHT: printf("Right"); break;
			case Impl::CornerTp::ACUTE: printf("Acute"); break;
			case Impl::CornerTp::REENTRANT: printf("Reentrant"); break;
			case Impl::CornerTp::ROUND: printf("Round"); break;
			case Impl::CornerTp::NEGLECTABLE: printf("Neglect."); break;
			default: printf("unknown");
		}
		std::cout<<std::endl;
	}
	std::cout<<"+++++++++++++++++++++++++++++++++++++"<<std::endl;
}


#endif
