#include "BsplineCurveEvaluator.h"
#include <assert.h>
#include "mat.h"
#include "vec.h"
#include "modelerapp.h"

#define NUM_SEG 20

void  BsplineCurveEvaluator::evaluateCurve(const std::vector<Point>& ptvCtrlPts,
	std::vector<Point>& ptvEvaluatedCurvePts,
	const float& fAniLength,
	const bool& bWrap) const
{
	ptvEvaluatedCurvePts.clear();
	const Mat4d convert = Mat4d(
		1, 4, 1, 0,
		0, 4, 2, 0,
		0, 2, 4, 0,
		0, 1, 4, 1) / 6.0;
	const Mat4d ker(
		-1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 3, 0, 0,
		1, 0, 0, 0);
	std::vector<Point> ptvCtrlPtsCpy;
	if (bWrap)
	{
		ptvCtrlPtsCpy.push_back(Point((ptvCtrlPts.end() - 2)->x - fAniLength,(ptvCtrlPts.end() - 2)->y));
		ptvCtrlPtsCpy.push_back(Point((ptvCtrlPts.end() - 1)->x - fAniLength, (ptvCtrlPts.end() - 1)->y));
		ptvCtrlPtsCpy.insert(ptvCtrlPtsCpy.end(), ptvCtrlPts.begin(), ptvCtrlPts.end());
		ptvCtrlPtsCpy.push_back(Point((ptvCtrlPts.begin())->x + fAniLength,(ptvCtrlPts.begin())->y));
		ptvCtrlPtsCpy.push_back(Point((ptvCtrlPts.begin() + 1)->x + fAniLength,(ptvCtrlPts.begin() + 1)->y));
	}
	else
	{
		ptvCtrlPtsCpy.assign(ptvCtrlPts.begin(), ptvCtrlPts.end());
		ptvCtrlPtsCpy.push_back(Point(ptvCtrlPts.front().x+fAniLength,ptvCtrlPts.front().y));
		ptvCtrlPtsCpy.insert(ptvCtrlPtsCpy.begin(), Point(ptvCtrlPts.back().x - fAniLength, ptvCtrlPts.back().y));
	}

	for (int iter = 0; iter < ((int)ptvCtrlPtsCpy.size() - 3); ++iter)
	{
		Vec4d ctrl_x(ptvCtrlPtsCpy[iter].x,ptvCtrlPtsCpy[iter + 1].x,
			ptvCtrlPtsCpy[iter + 2].x,ptvCtrlPtsCpy[iter + 3].x);
		Vec4d ctrl_y(ptvCtrlPtsCpy[iter].y,ptvCtrlPtsCpy[iter + 1].y,
			ptvCtrlPtsCpy[iter + 2].y,ptvCtrlPtsCpy[iter + 3].y);
		Vec4d convert_x = convert * ctrl_x;
		Vec4d convert_y = convert * ctrl_y;
		std::vector<Point> ptvEvalPts;
		for (int i = 0; i < NUM_SEG; ++i)
		{
			double t = i / (double)NUM_SEG;
			Vec4d para(t*t*t, t*t, t, 1);
			Point eval_pt((float)(para*ker*convert_x), (float)(para*ker*convert_y));
			ptvEvalPts.push_back(eval_pt);
		}
		ptvEvaluatedCurvePts.insert(ptvEvaluatedCurvePts.end(), ptvEvalPts.begin(), ptvEvalPts.end());
	}
	if (!bWrap)
	{
		ptvEvaluatedCurvePts.push_back(Point(0, ptvCtrlPts.front().y));
		ptvEvaluatedCurvePts.push_back(Point(fAniLength, ptvCtrlPts.back().y));
	}
}