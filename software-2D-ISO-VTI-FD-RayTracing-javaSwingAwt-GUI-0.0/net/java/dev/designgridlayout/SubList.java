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

import java.util.AbstractList;
import java.util.List;
import java.util.RandomAccess;

// Simplified sublist implementation tailored for managing RowItem lists in
// SubGrid/GridRow.
// This implementation is definitly not complete for external use but perfectly
// fulfills its goals in the required context.
final class SubList extends AbstractList<RowItem> implements RandomAccess
{
	SubList(List<RowItem> source)
	{
		_source = source;
		_from = source.size();
	}
	
	@Override public void add(int index, RowItem element)
	{
		_source.add(_from + index, element);
		_size++;
	}

	@Override public RowItem get(int index)
	{
		return _source.get(_from + index);
	}

	@Override public int size()
	{
		return _size;
	}

	private final List<RowItem> _source;
	private final int _from;
	private int _size = 0;
}
