//  Copyright 2005-2011 Jason Aaron Osgood, Jean-Francois Poilpret
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package net.java.dev.designgridlayout;

import java.awt.Container;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.LayoutStyle.ComponentPlacement;

final class ComponentHelper
{
	private ComponentHelper()
	{
	}
	
	static int maxValues(Iterable<? extends IItem> items, IExtractor extractor)
	{
		int max = 0;
		for (IItem item: items)
		{
			max = Math.max(max, extractor.value(item));
		}
		return max;
	}
	
	static int sumValues(Iterable<? extends IItem> items, IExtractor extractor)
	{
		int sum = 0;
		for (IItem item: items)
		{
			sum += extractor.value(item);
		}
		return sum;
	}
	
	static int unrelhgap(List<? extends IItem> items, Container parent)
	{
		return hgap(
			ComponentGapsHelper.instance(), items, parent, ComponentPlacement.UNRELATED);
	}
	
	static int hgap(List<? extends IItem> items, Container parent)
	{
		return hgap(
			ComponentGapsHelper.instance(), items, parent, ComponentPlacement.RELATED);
	}
	
	static int hgap(JComponent first, List<? extends IItem> items, Container parent)
	{
		ComponentGapsHelper helper = ComponentGapsHelper.instance();
		int hgap = 0;
		if (first != null && !items.isEmpty())
		{
			hgap = hgap(helper, first, items.get(0).component(), parent, 
				ComponentPlacement.RELATED);
		}
		return Math.max(hgap, hgap(helper, items, parent, ComponentPlacement.RELATED));
	}

	static boolean isFixedHeight(
		HeightGrowPolicy policy, Iterable<? extends IRowItem> items)
	{
		for (IRowItem item: items)
		{
			if (policy.canGrowHeight(item.component()) && item.isLastSpanRow())
			{
				return false;
			}
		}
		return true;
	}

	static private int hgap(ComponentGapsHelper helper, List<? extends IItem> items, 
		Container parent, ComponentPlacement placement)
	{
		int hgap = 0;
		int size = items.size() - 1;
		for (int nth = 0; nth < size; nth++)
		{
			JComponent left = items.get(nth).component();
			JComponent right = items.get(nth + 1).component();
			int gap = hgap(helper, left, right, parent, placement);
			hgap = Math.max(hgap, gap);
		}
		return hgap;
	}
	
	static private int hgap(ComponentGapsHelper helper, JComponent left, JComponent right, 
		Container parent, ComponentPlacement placement)
	{
		return helper.getHorizontalGap(left, right, placement, parent);
	}
}
