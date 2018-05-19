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

import java.util.ArrayList;
import java.util.List;

import javax.swing.JComponent;

abstract class AbstractNonGridRow extends AbstractRow implements INonGridRow
{
	@Override public INonGridRow add(JComponent... children)
	{
		for (JComponent component: children)
		{
			checkAddedComponent(component);
			_items.add(new NonGridRowItem(component));
			parent().add(component);
		}
		return this;
	}

	@Override public INonGridRow addMulti(JComponent... children)
	{
		for (JComponent component: children)
		{
			checkAddedComponent(component);
		}
		return add(Componentizer.create().fixedPref(children).component());
	}

	@Override public INonGridRow indent()
	{
		return indent(1);
	}
	
	@Override public INonGridRow indent(int n)
	{
		if (n >= 0)
		{
			_indent = ComponentGapsHelper.instance().getHorizontalIndent() * n;
		}
		return this;
	}
	
	@Override public INonGridRow fill()
	{
		_fill = true;
		return this;
	}

	@Override public INonGridRow withOwnRowWidth()
	{
		_ownRowWidth = true;
		return this;
	}

	@Override List<NonGridRowItem> items()
	{
		return _items;
	}

	@Override int totalNonGridWidth(int hgap, int unrelhgap)
	{
		int count = _items.size();
		int totalWidth = _compWidth * count + (hgap * (count - 1)) + _indent;
		return totalWidth;
	}

	@Override int componentNonGridWidth()
	{
		return (_ownRowWidth ? 0 : actualComponentNonGridWidth());
	}
	
	private int actualComponentNonGridWidth()
	{
		return ComponentHelper.maxValues(_items, PrefWidthExtractor.INSTANCE);
	}

	@Override void forceComponentNonGridWidth(int width)
	{
		_compWidth = (width > 0 && !_ownRowWidth ? width : actualComponentNonGridWidth());
	}

	//CSOFF: ParameterNumber
	@Override int layoutRow(LayoutHelper helper, int left, int hgap, int gridgap, 
		int unrelhgap, int rowWidth, int gridsWidth, List<Integer> labelsWidth)
	{
		// Calculate various needed widths & origin
		int x = left + _indent;
		int count = _items.size();
		int width = _compWidth;

		int leftFiller = width;
		int rightFiller = width;
		if (!_fill)
		{
			int usedWidth = width * count + ((count - 1) * hgap) + _indent;
			x += xOffset(rowWidth, usedWidth);
		}
		else
		{
			int availableWidth = rowWidth - _indent - (count - 1) * hgap;
			leftFiller = leftFiller(count, width, availableWidth);
			rightFiller = rightFiller(count, width, availableWidth);
		}

		return layoutRow(helper, x, hgap, width, leftFiller, rightFiller);
	}
	//CSON: ParameterNumber

	protected int layoutRow(LayoutHelper helper, int left, int hgap, int width, 
		int leftFiller, int rightFiller)
	{
		int x = left;
		int count = _items.size();
		int i = 0;
		int actualHeight = 0;
		for (NonGridRowItem item: _items)
		{
			int compWidth;
			if (i == 0)
			{
				compWidth = leftFiller;
			}
			else if (i == count - 1)
			{
				compWidth = rightFiller;
			}
			else
			{
				compWidth = width;
			}
			actualHeight = Math.max(actualHeight, helper.setSizeLocation(
				item.component(), x, compWidth, height(), baseline()));
			x += compWidth + hgap;
			i++;
		}
		return actualHeight;
	}

	abstract protected int xOffset(int rowWidth, int usedWidth);

	abstract protected int leftFiller(int count, int width, int availableWidth);

	abstract protected int rightFiller(int count, int width, int availableWidth);

	private final List<NonGridRowItem> _items = new ArrayList<NonGridRowItem>();
	private boolean _fill = false;
	private boolean _ownRowWidth = false;
	private int _compWidth = 0;
	private int _indent = 0;
}
