#include "CatmullRomCurveEvaluator.h"
#include <assert.h>
#include "mat.h"
#include "vec.h"

#define NUM_SEG 30

void CatmullRomCurveEvaluator::evaluateCurve(const std::vector<Point>& ptvCtrlPts,
	std::vector<Point>& ptvEvaluatedCurvePts,
	const float& fAniLength,
	const bool& bWrap) const
{
	ptvEvaluatedCurvePts.clear();
	const Mat4d ker = Mat4d(
		-1, 3, -3, 1,
		2, -5, 4, -1,
		-1, 0, 1, 0,
		0, 2, 0, 0) / 2.0;
	std::vector<Point> ptvCtrlPtsCpy;
	if (bWrap)
	{
		ptvCtrlPtsCpy.push_back(Point((ptvCtrlPts.end() - 2)->x - fAniLength, (ptvCtrlPts.end() - 2)->y));
		ptvCtrlPtsCpy.push_back(Point((ptvCtrlPts.end() - 1)->x - fAniLength, (ptvCtrlPts.end() - 1)->y));
		ptvCtrlPtsCpy.insert(ptvCtrlPtsCpy.end(), ptvCtrlPts.begin(), ptvCtrlPts.end());
		ptvCtrlPtsCpy.push_back(Point((ptvCtrlPts.begin())->x + fAniLength, (ptvCtrlPts.begin())->y));
		ptvCtrlPtsCpy.push_back(Point((ptvCtrlPts.begin() + 1)->x + fAniLength, (ptvCtrlPts.begin() + 1)->y));
	}
	else
	{
		ptvCtrlPtsCpy.assign(ptvCtrlPts.begin(), ptvCtrlPts.end());
		ptvCtrlPtsCpy.push_back(Point(ptvCtrlPts.front().x + fAniLength, ptvCtrlPts.front().y));
		ptvCtrlPtsCpy.insert(ptvCtrlPtsCpy.begin(), Point(ptvCtrlPts.back().x - fAniLength, ptvCtrlPts.back().y));
	}

	for (int iter = 0; iter < ((int)ptvCtrlPtsCpy.size() - 3); ++iter)
	{
		Vec4d ctrl_x = Vec4d(ptvCtrlPtsCpy[iter].x,ptvCtrlPtsCpy[iter + 1].x,
			ptvCtrlPtsCpy[iter + 2].x,ptvCtrlPtsCpy[iter + 3].x);
		Vec4d ctrl_y = Vec4d(ptvCtrlPtsCpy[iter].y,ptvCtrlPtsCpy[iter + 1].y,
			ptvCtrlPtsCpy[iter + 2].y,ptvCtrlPtsCpy[iter + 3].y);
		std::vector<Point> ptvEvalPts;
		for (int i = 1; i < NUM_SEG; ++i)
		{
			double t = i / (double)NUM_SEG;
			Vec4d para = Vec4d(t*t*t, t*t, t, 1);
			Point eval_pt = Point((float)(para * ker * ctrl_x), (float)(para * ker * ctrl_y));
			if (eval_pt.x > ptvCtrlPtsCpy[iter + 1].x && eval_pt.x < ptvCtrlPtsCpy[iter + 2].x
				&& (ptvEvalPts.empty() || eval_pt.x > ptvEvalPts.back().x))
				ptvEvalPts.push_back(eval_pt);
		}
		ptvEvaluatedCurvePts.push_back(ptvCtrlPtsCpy[iter + 1]);
		ptvEvaluatedCurvePts.insert(ptvEvaluatedCurvePts.end(), ptvEvalPts.begin(), ptvEvalPts.end());
		ptvEvaluatedCurvePts.push_back(ptvCtrlPtsCpy[iter + 2]);
	}

	if (!bWrap)
	{
		ptvEvaluatedCurvePts.push_back(Point(0, ptvCtrlPts.front().y));
		ptvEvaluatedCurvePts.push_back(Point(fAniLength, ptvCtrlPts.back().y));
	}
}